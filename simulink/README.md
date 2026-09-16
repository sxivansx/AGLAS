# Simulink model: adaptive gust load alleviation

Real-time simulation of the AGLAS wing with an LQG baseline and an MRAC
adaptive augmentation, in state-space form.

## Run it

```matlab
cd simulink
aglas_sim_setup          % design the controller, put matrices in the base workspace
build_gla_model          % generate and open aglas_gla.slx
compare_simulink         % runs the model and checks it against the .m reference
```

`compare_simulink` runs the simulation itself, so there is no need to call
`sim` separately. It is the important step: run it before trusting any result.

## Logged signals

All six are prefixed `aglas_`: `aglas_moment`, `aglas_tip`, `aglas_delta_cmd`,
`aglas_theta`, `aglas_u_ad`, `aglas_gust`.

The prefix is not cosmetic. An earlier version logged to a variable called
`moment`, and because recent MATLAB returns a `Simulink.SimulationOutput`
object rather than writing To Workspace variables into the base workspace, the
retrieval found no variable and resolved to MATLAB's built-in `moment()`
function instead. The error, `Moment requires two inputs`, said nothing about
the actual problem. Retrieval now reads the `SimulationOutput` object first and
falls back to the base workspace, and reports the missing signal by name if
neither has it.

## What the model contains

| Block | Contents |
|---|---|
| Plant | 8-state aeroservoelastic model: 3 mass-normalised modes plus a 2nd-order actuator. Inputs are surface deflection and gust velocity. |
| Estimator | 9-state Kalman filter. Eight structural states plus a gust state, so it can tell a gust apart from the wing's response to one. |
| Kgain | LQR gain. Acting on the gust estimate too, which is what makes it disturbance feedforward rather than a purely reactive loop. |
| RefModel | The nominal closed loop. MRAC drives the aircraft towards this. |
| AdaptLaw | Parameter update with normalisation, sigma-modification and a soft norm barrier. |
| RateLim, Sat | Actuator limits, 150 deg/s and 20 deg. |

Solver is fixed-step `ode4` at 2e-4 s. The rate limiter is a stateful
nonlinearity whose output depends on the previously accepted step, so a
variable-step solver that rejects and retries steps can corrupt it.

## Switching between LQG and MRAC

```matlab
AGLAS.adapt_on = 0;   % plain LQG
AGLAS.adapt_on = 1;   % MRAC augmentation
```

No rebuild needed. It is a gain on the adaptive branch.

## Flying off-design

The fourth argument to `aglas_sim_setup` scales the *true* control
effectiveness while leaving the controller's model untouched, which is how you
simulate a degraded, iced or partly jammed surface:

```matlab
aglas_sim_setup('baseline', 0.3, 5, 0.4);   % surface is 40 % as effective
sim('aglas_gla')
```

The third argument is the MRAC adaptation rate. It is sensitive: see the note
on results below.

## Expected results at the nominal condition

From the pure-MATLAB reference, which the model should reproduce:

| Case | Peak root moment | Reduction | Peak deflection |
|---|---|---|---|
| Open loop | 145.2 kN m | — | 0 deg |
| LQG | 68.4 kN m | 52.9 % | 17.9 deg |
| MRAC, gamma = 0.5 | 71.4 kN m | 50.8 % | 17.3 deg |

MRAC being slightly *worse* than LQG at the design point is correct, not a
fault. There is no model error to adapt to there, so the adaptive term can only
add activity. The honest reading of these numbers is in the top-level control
documentation: on this aircraft, across a realistic 35 to 200 m/s envelope,
fixed-gain LQG stays stable and adaptation does not earn its complexity.

## Interpreting the comparison

`compare_simulink` reports four numbers. The ones that matter are the peak
differences and the waveform RMS. A large worst-pointwise difference alongside
small peak differences is normal: it means the two traces agree on the load the
structure is sized by, and differ slightly in timing during a lightly damped
secondary oscillation.

One difference between the two is structural and permanent. The `.m` simulation
computes the control once per step and holds it across the four RK4 stages,
which is what a digital controller sampling at 5 kHz actually does. Simulink
integrates the controller as part of one continuous system, so its command
varies within the step.

If the peaks disagree, set `AGLAS.adapt_on = 0` and re-run. That removes both
adaptive integrators. If plain LQG then agrees tightly, the problem is in the
adaptive path, not in the plant or estimator wiring.

## Verification status

The `.m` path is covered by the project test suite and was used to produce
every number above. **The Simulink model itself has not been executed**; it was
generated against the documented block API but no Simulink licence was
available where it was written. `compare_simulink` exists precisely for that
reason. A block diagram is easy to get subtly wrong, a swapped Mux order, a
Selector off by one, a sign on a Sum, and none of those raise an error. They
just give a plausible wrong answer.

If `compare_simulink` reports a disagreement, check in this order:

1. `AdaptIn` Mux order: it must be `[xhat(9); xm(8); theta(8)]`.
2. `UadIn` Mux order: `[theta(8); xhat_struct(8)]`.
3. `PickStruct` indices `1:8`, `PickGust` index `9`.
4. The `SumU` signs, both `+`.
5. That `Kgain` has `Multiplication` set to `Matrix(K*u)`.
