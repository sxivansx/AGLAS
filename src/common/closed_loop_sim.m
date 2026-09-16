function out = closed_loop_sim(cfg, P, ctrl, t, wg)
%CLOSED_LOOP_SIM  Time simulation of the wing with a gust and a controller.
%
%   out = CLOSED_LOOP_SIM(cfg, P, ctrl, t, wg)
%
%   ctrl.type selects the law:
%     'none'  open loop, the surface stays at zero
%     'lqr'   full state feedback, u = -K*x. Not implementable, since modal
%             states are not measurable, but it bounds what any output-feedback
%             design can achieve and so is the right yardstick for the others.
%     'lqg'   the implementable law: Kalman estimator driving the same gain
%     'pid'   classical single loop on root strain, for comparison
%
%   Fixed-step RK4 rather than ode45. The actuator rate limiter is a stateful
%   nonlinearity whose output depends on the previous accepted command, so it
%   is not a function of the current state alone. An adaptive solver that
%   rejects and retries steps would feed it a non-monotonic time sequence and
%   corrupt it. A fixed step keeps the limiter honest.
%
%   Actuator limits are applied to the command: amplitude clipped at
%   cfg.control.delta_max and slew clipped at cfg.control.rate_max, ahead of
%   the second-order actuator already present in the plant. Peak actual
%   deflection and rate are reported so it is visible whether the physical
%   surface was driven past its limits.
%
%   Output fields
%     t, x, delta_cmd, delta, delta_rate
%     moment, tip, strain, accel      performance and measurement histories
%     peak                            struct of peak magnitudes
%     saturated                       fraction of samples on a limit

    t  = t(:).';
    wg = wg(:).';
    n  = numel(t);
    dt = t(2) - t(1);
    if max(abs(diff(t) - dt)) > 1e-12
        error('closed_loop_sim:nonUniform', 'Time vector must be uniformly spaced.');
    end

    ns = P.n_states;
    x  = zeros(ns, n);
    u  = zeros(1, n);

    % The estimator carries an extra gust state, so its model is larger than
    % the plant. 'lqr' uses the true state augmented with the true gust, which
    % is the ideal every implementable law is measured against.
    if isfield(ctrl, 'n_aug') && ~isempty(ctrl.n_aug)
        na = ctrl.n_aug;
    else
        na = ns;
    end
    xh = zeros(na, n);

    % MRAC extra states: reference model and adaptive parameters.
    is_mrac = strcmpi(ctrl.type, 'mrac');
    if is_mrac
        xm    = zeros(ns, n);
        theta = zeros(ns, n);
    else
        xm = []; theta = [];
    end

    d_max = deg2rad(cfg.control.delta_max);
    r_max = deg2rad(cfg.control.rate_max);

    use_limits = true;
    if isfield(ctrl, 'ideal_actuator') && ctrl.ideal_actuator
        use_limits = false;
    end

    % --- explicit RK4 stability guard --------------------------------------
    % Classical RK4 is stable only while |lambda|*dt stays below about 2.78. A
    % Kalman filter tuned with an over-confident sensor model can easily have
    % poles hundreds of times faster than the structure, and the simulation
    % then diverges silently while still producing plausible-looking peak
    % numbers. An earlier run of this file did exactly that: it reported a 97 %
    % load reduction that was really an estimator overflowing to 1e300 and
    % slamming the surface between its stops. Fail loudly instead.
    lam_plant = max(abs(eig(P.A)));
    lam_max   = lam_plant;
    if isfield(ctrl, 'Ae') && ~isempty(ctrl.Ae) && ...
       (strcmpi(ctrl.type, 'lqg') || strcmpi(ctrl.type, 'mrac'))
        lam_max = max(lam_max, max(abs(eig(ctrl.Ae - ctrl.L*ctrl.Ce))));
    end
    if lam_max * dt > 2.5
        error('closed_loop_sim:stepTooLarge', ...
            ['Time step %.2e s cannot integrate a loop whose fastest pole is ' ...
             '%.1f rad/s (%.1f Hz): RK4 needs |lambda|*dt below about 2.78, ' ...
             'and this gives %.2f. Reduce dt below %.2e s, or slow the ' ...
             'estimator by lowering cfg.sensor.gust_psd.'], ...
             dt, lam_max, lam_max/(2*pi), lam_max*dt, 2.5/lam_max);
    end

    n_sat_pos  = 0;
    n_sat_rate = 0;
    u_prev = 0;
    bound  = 1e6;

    % PID state
    e_int = 0;
    e_prev = 0;

    for k = 1:n-1
        xk  = x(:,k);
        xhk = xh(:,k);
        wk  = wg(k);

        % ---- measurement (used by output-feedback laws) -------------------
        yk = P.C*xk + P.D*u_prev + P.Dw*wk;

        % ---- control law ---------------------------------------------------
        switch lower(ctrl.type)
            case 'none'
                uc = 0;
            case 'lqr'
                % Full state plus perfect knowledge of the gust.
                uc = -ctrl.K * [xk; wk];
            case 'lqg'
                uc = -ctrl.K * xhk;
            case 'mrac'
                % Baseline LQG command plus the adaptive increment.
                phi = xhk(1:ns);
                uc  = -ctrl.K * xhk + theta(:,k).' * phi;
            case 'pid'
                % Single loop on root strain. Sign chosen so a positive strain,
                % meaning the wing is bending up, commands the surface to
                % reduce lift.
                e = yk(1);
                e_int = e_int + e*dt;
                de = (e - e_prev)/dt;
                e_prev = e;
                uc = -(ctrl.kp*e + ctrl.ki*e_int + ctrl.kd*de);
            otherwise
                error('closed_loop_sim:unknownType', ...
                      'Unknown controller type "%s".', ctrl.type);
        end

        % ---- actuator command limits ---------------------------------------
        % Position limit first, then slew limit measured from the position the
        % surface actually reached. That ordering is what a real actuator does:
        % it cannot be commanded past its stops, and from wherever it currently
        % sits it cannot move faster than its hydraulic rate.
        %
        % The reverse order, slew then clip, is subtly different and the two
        % disagree exactly while the surface is on a stop. A rate limiter
        % holding its own pre-saturation output as state keeps winding past the
        % stop and then has to unwind before the surface moves back, which is
        % actuator windup. Doing it in this order removes the windup and, just
        % as usefully, matches the Simulink model block for block.
        if use_limits
            if abs(uc) > d_max
                uc = sign(uc)*d_max;
                n_sat_pos = n_sat_pos + 1;
            end
            du = uc - u_prev;
            if abs(du) > r_max*dt
                uc = u_prev + sign(du)*r_max*dt;
                n_sat_rate = n_sat_rate + 1;
            end
        end
        u(k) = uc;
        u_prev = uc;

        % ---- RK4 on the plant ----------------------------------------------
        wk2 = interp_w(t, wg, t(k) + dt/2);
        wk3 = wg(k+1);
        f = @(xx, ww) P.A*xx + P.B*uc + P.Bw*ww;
        k1 = f(xk,             wk);
        k2 = f(xk + dt/2*k1,   wk2);
        k3 = f(xk + dt/2*k2,   wk2);
        k4 = f(xk + dt*k3,     wk3);
        x(:,k+1) = xk + dt/6*(k1 + 2*k2 + 2*k3 + k4);

        if ~all(isfinite(x(:,k+1))) || max(abs(x(:,k+1))) > bound
            error('closed_loop_sim:diverged', ...
                ['Plant states diverged at t = %.4f s (max |x| = %.3e). The ' ...
                 'controller is destabilising the plant, or the step is too ' ...
                 'large.'], t(k+1), max(abs(x(:,k+1))));
        end

        % ---- reference model and adaptive parameter update ------------------
        if is_mrac
            w_hat = xhk(ctrl.idx_gust);
            e     = xhk(1:ns) - xm(:,k);

            % RK4 on the reference model, to match the ode4 solver the
            % Simulink model uses. Forward Euler here was a visible source of
            % disagreement between the two.
            fm = @(xx) ctrl.Am*xx + P.Bw*w_hat;
            r1 = fm(xm(:,k));
            r2 = fm(xm(:,k) + dt/2*r1);
            r3 = fm(xm(:,k) + dt/2*r2);
            r4 = fm(xm(:,k) + dt*r3);
            xm(:,k+1) = xm(:,k) + dt/6*(r1 + 2*r2 + 2*r3 + r4);

            % theta_hat_dot = -Gamma*( Phi*(B'*P_lyap*e)/(1+Phi'*Phi) + sigma*theta )
            phi_a = xhk(1:ns);
            s_e   = ctrl.PB.' * e;                    % scalar
            if ctrl.normalise
                s_e = s_e / (1 + phi_a.'*phi_a);
            end
            dth   = -ctrl.Gamma * (phi_a*s_e + ctrl.sigma*theta(:,k));
            th_n  = theta(:,k) + dt*dth;

            % Hard projection keeps the parameters bounded even if the leakage
            % term is out-tuned by a persistent unmatched disturbance.
            nrm = norm(th_n);
            if nrm > ctrl.theta_max
                th_n = th_n * (ctrl.theta_max / nrm);
            end
            theta(:,k+1) = th_n;
        end

        % ---- RK4 on the estimator (it does not know the gust) ---------------
        if strcmpi(ctrl.type, 'lqg') || is_mrac
            g = @(xx, yy) ctrl.Ae*xx + ctrl.Be*uc + ctrl.L*(yy - ctrl.Ce*xx - ctrl.De*uc);
            y_mid = P.C*(xk + dt/2*k1) + P.D*uc + P.Dw*wk2;
            y_end = P.C*x(:,k+1)       + P.D*uc + P.Dw*wk3;
            m1 = g(xhk,            yk);
            m2 = g(xhk + dt/2*m1,  y_mid);
            m3 = g(xhk + dt/2*m2,  y_mid);
            m4 = g(xhk + dt*m3,    y_end);
            xh(:,k+1) = xhk + dt/6*(m1 + 2*m2 + 2*m3 + m4);
            if ~all(isfinite(xh(:,k+1))) || max(abs(xh(:,k+1))) > bound
                error('closed_loop_sim:estimatorDiverged', ...
                    ['Estimator diverged at t = %.4f s (max |xhat| = %.3e). ' ...
                     'The filter is too fast for this step size.'], ...
                     t(k+1), max(abs(xh(:,k+1))));
            end
        end
    end
    u(n) = u(n-1);

    % ------------------------------------------------------------- outputs
    z      = P.Cz * x;
    meas   = P.C  * x + P.D * u + P.Dw * wg;

    out.t          = t;
    out.x          = x;
    out.xhat       = xh;
    out.xm         = xm;
    out.theta      = theta;
    out.delta_cmd  = u;
    out.delta      = x(P.idx.delta, :);
    out.delta_rate = x(P.idx.ddot,  :);
    out.moment     = z(1,:);
    out.tip        = z(2,:);
    out.strain     = meas(1,:);
    out.accel      = meas(2,:);
    out.wg         = wg;

    out.peak.moment      = max(abs(out.moment));
    out.peak.tip         = max(abs(out.tip));
    out.peak.delta_deg   = rad2deg(max(abs(out.delta)));
    out.peak.rate_deg_s  = rad2deg(max(abs(out.delta_rate)));
    out.peak.cmd_deg     = rad2deg(max(abs(out.delta_cmd)));

    out.saturated.position = n_sat_pos  / max(1, n-1);
    out.saturated.rate     = n_sat_rate / max(1, n-1);
    out.limits.delta_deg   = cfg.control.delta_max;
    out.limits.rate_deg_s  = cfg.control.rate_max;
    out.type = ctrl.type;
end

% ------------------------------------------------------------------------
function w = interp_w(t, wg, tq)
    if tq <= t(1)
        w = wg(1);
    elseif tq >= t(end)
        w = wg(end);
    else
        dt = t(2) - t(1);
        i0 = min(max(floor((tq - t(1))/dt) + 1, 1), numel(t)-1);
        a  = (tq - t(i0))/dt;
        w  = (1-a)*wg(i0) + a*wg(i0+1);
    end
end
