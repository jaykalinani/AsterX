# EOS repair series

Review snapshot: 2026-09-29. Repository: `repos/AsterX`; branch: `eos_repairs`.
The local `dev` and existing `origin/dev` references both remain at
`d91110629c68ee2942ca916a12e131438b70c0b6`. No remote reference was fetched or
updated for this review. No production parfile or executable was changed by
this series.

Status: source changes and opt-in tests are committed for review. **Nothing
has been compiled or rebuilt, and the new C++ tests have not been run.**
Static review is not evidence of successful compilation, GPU execution, or
improved stellar evolution. The user will perform the rebuild.

## Scope and numbering

The supported targets of this series are ideal gas and tabulated 3-D EOSs.
Hybrid EOS validation is excluded. Its dispatch remains in the code, but
shared-template changes must not be taken as a claim that hybrid is tested.

The proposed list had different numbering in its summary and detailed
sections. This report follows the detailed list: 13 is atmosphere thresholds,
14 is diagnostics/tests, and 15 is the added EOS-agnostic C2P repair.
Positivity-preserving flux limiting remains a separate, deferred experiment.

| # | Repair and implementation | Main commits and review follow-ups |
|---|---|---|
| 1 | Physical tabulated energy range; shift only for interpolation; bound both inverse endpoints | `045b7579` |
| 2 | Common temperature-, energy-, pressure-, and enthalpy-state closure | `ca53dd48`, `5859a457`, `23a2e38f`, `4c92b59d` |
| 3 | Shared graded atmosphere builder; preserve selected polytropic/piecewise-polytropic cold EOS | `e9e2df2c`, `803cffe6`, `cf3f9686` |
| 4 | EOS-consistent initial data and primitive checks | `f0977e22`, `2442f501`, `cd2e7d70` |
| 5 | Complete atmosphere resets; complete optional neighbour repair | `803cffe6`, `62650754` |
| 6 | Separate EOS domain, atmosphere state, and cutoff policy without renaming established fields | `9a99a326`, `10df9d89`, `2895402d` |
| 7 | EOS-aware tau repair inside the active C2P path | `e964570a`, `f682f148` |
| 8 | Shared post-C2P closure, local energy bounds, bounded magnetic heating, finite-state rejection | `3f2c77d9`, `f267d607`, `cf8a7d39`, `25e0d49c`, `ac208ef6` |
| 9 | EOS-valid reconstructed face states with local graded atmosphere | `6bd319bb`, `cd2e7d70` |
| 10 | Evolved kappa consistently used in repaired states, face fluxes, and beta-floor changes | `f0977e22`, `f267d607`, `6bd319bb`, `9a649102` |
| 11 | Sound speed from the same face state; degenerate HLLE guard | `6bd319bb`, `624f2fa1`, `96b57af7`, `791f4692` |
| 12 | Startup bounds/mode validation and effective-parameter reporting | `ed450532`, `2895402d`, `cf3f9686`, `cd2e7d70` |
| 13 | Explicit opt-in matching of face and cell atmosphere thresholds | `801f2b0d`, `b54b7b43` |
| 14 | Optional EOS-call/event counters and thorn-local test suites | `ca24bdab`, `96b57af7`, `426f08af`, `926868db`, `224787c2`, `a5c529bb`, `791f4692`, `55670f80`, `8f34d4be` |
| 15 | Remove GammaIdealFluid from C2Ps; general-EOS Noble Jacobian; RePrimAnd negative-energy enthalpy bound; remove standalone RPA file | `25005ffe`, `a24f4b9a`, `3adada4e`, `53d6d61b`, `2ee989aa`, `c84348f3` |

The earlier printf correction is retained in `e796ebff`: the point index is
printed as three integer coordinates, rather than passing a vector to `%d`.
The intermediate cutoff rename was explicitly undone by `10df9d89`.
`atmo.rho_cut`, `rho_abs_min`, `tauFluid_atmo`, and existing solver variable
names remain. The requested EOS diagnostic parameter is exactly
`EOSX::eos_call_diagnostics_every`, not `EOSX::eos_call_every`.

Use `git log --reverse --oneline dev..eos_repairs` for the complete ordered
commit list and `git diff dev..eos_repairs -- <thorn>` for a thorn-level review.
The deleted `AsterX/src/con2prim_rpa.cxx` is recoverable from `dev` or the parent
of `c84348f3`; the active `Con2PrimFactory` RePrimAnd solver remains.

## One thermodynamic authority

Implementation: [EOSX/src/thermo_state.hxx](EOSX/src/thermo_state.hxx).

| Mode | Independent inputs | Dependent outputs | Intended use |
|---|---|---|---|
| Temperature | rho, T, Ye | eps, P, kappa, cs2 | Tabulated ID, temperature reconstruction/floors; also supported for ideal gas |
| Energy | rho, eps, Ye | T, P, kappa, cs2 | Recovered states, cold atmosphere matching, saved C2P seeds |
| Pressure | rho, P, Ye | eps, T, kappa, cs2 | Ideal gas, where inversion is analytic |
| Enthalpy | rho, h, Ye | eps, P, T, kappa, derivatives | General-EOS Noble iteration |

The temperature path bounds rho/T/Ye and then derives the rest. It does not
independently impose a positive eps floor. The energy path bounds rho/Ye,
obtains the energy interval at that rho/Ye, bounds eps, then derives the rest.
The pressure path first obtains energy and uses the energy closure, so a
pressure outside the realizable EOS domain need not be returned unchanged.

These helpers bound inputs; they are not a general input-error reporting
interface. In particular, the fmin/fmax bound helper alone does not diagnose
NaNs. Recovery finalization and face/primitive checks explicitly reject or
reset non-finite independent states. A malformed table is not made valid by
clamping its input coordinates.

For ideal gas, Gamma remains inside the EOS implementation, where it belongs:

```
P = (Gamma - 1) rho eps
T = (Gamma - 1) particle_mass eps
kappa = P / rho^Gamma
cs2 = Gamma (Gamma - 1) eps / (1 + Gamma eps)
```

The C2P algorithms no longer read a GammaIdealFluid member or an EOS gamma
field. The pressure derivative and inversion capabilities come from the EOS.

## Bounds, atmosphere, and reset policies

| Quantity | Owner/source | Where it is used | Important distinction |
|---|---|---|---|
| `rgrho`, `rgtemp`, `rgye` | Table axes for Tabulated3d; configured/derived bounds for ideal gas | Closure, reconstruction, recovery, atmosphere construction | EOS validity, not the atmosphere profile |
| `rgeps` | Physical global table energy extrema; ideal-gas configured range, with its existing causality cap | Global checks, table inversion guard, conservative admissibility trigger | Not the energy interval at every rho/Ye |
| `range_eps_from_rho_ye` | Table temperature endpoints at fixed rho/Ye; ideal gas returns its configured range | Recovery, energy closure, repaired tau target | Can have a negative minimum for the table |
| `rho_abs_min`, `r_atmo`, `n_rho_atmo` | Atmosphere policy | Local atmosphere density before EOS-domain bounding | Not aliases for the table's minimum density |
| `t_atmo`, `n_temp_atmo` | Temperature-primary atmosphere | Graded T, then EOS closure | Ignored when constructing a pressure-primary or cold-matched atmosphere |
| `p_atmo`, `n_press_atmo` | Pressure-primary atmosphere | Graded P, then supported EOS inversion | Not a generic tabulated pressure inversion |
| `Ye_atmo` | Atmosphere composition | Bounded by `rgye`, then used for the complete atmosphere | Not the minimum allowed Ye of evolved cells |
| `atmo.eps_atmo`, `press_atmo`, `temp_atmo`, `entropy_atmo` | Shared atmosphere builder | Exact resets and selected local thermal-floor policy | Derived together; eps_atmo is not a universal C2P energy floor |
| `eps_atmo` parameter | Legacy accepted input | Startup explains that the builder derives energy | No independent energy constraint is added |
| `atmo.rho_cut` | `atmo.rho_atmo * (1 + atmo_tol)` | Cell atmosphere classification | A reset threshold, not an EOS bound |
| `ReconX::recon_thresh` | Existing face policy | Local face-side atmosphere density multiplier | Default 0 remains unchanged |
| `ReconX::recon_use_atmo_tol` | New opt-in switch, default no | Uses `1 + atmo_tol` for face threshold | Does not change existing parfiles silently |
| `tauFluid_atmo` | Numerical conservative-repair margin | Only after the global tau admissibility test fails | Undensitized energy density, not specific energy or atmosphere thermodynamics |

Implementation: [Con2PrimFactory/src/atmo.hxx](Con2PrimFactory/src/atmo.hxx).
For each requested spatial radius, density remains

```
rho_requested = rho_abs_min                       (r <= r_atmo)
rho_requested = rho_abs_min * (r_atmo/r)^n_rho_atmo (r >  r_atmo)
rho_atm = clamp(rho_requested, EOS rho range)
```

Temperature and pressure have their analogous existing grading laws in the
modes where they are authoritative. The builder bounds the independent state
and computes all dependent thermodynamics. At flux interfaces, it is evaluated
at the **two original adjacent cell-center positions**, separately, not once
at the face midpoint. This preserves the pre-existing sampling geometry.

Cold matching takes eps from the selected cold EOS at the graded density,
then closes with the evolution EOS. For `poly_gamma = gl_gamma = 2`, away
from EOS-domain saturation, the expected scalings are rho/eps/T proportional
to r^-6 and P proportional to r^-12 when `n_rho_atmo=6`. For unequal Gammas,
energy matching does not imply pressure matching; startup reports the mismatch.

The builder now also serves ideal-gas face atmospheres. Therefore preserving
the intended cold grading is not a claim of bitwise equivalence to the old
face path, which could use different thermal inputs from cell recovery.

## Physical energy versus interpolation shift

Specific internal energy is measured relative to the EOS's rest-mass
convention. Nuclear binding can make tabulated physical eps negative without
making temperature, pressure, or total rest-frame energy negative. The
GRMHD formulas here contain `rho * (1 + eps)` and `h = 1 + eps + P/rho`.
Replacing a valid negative eps with zero adds physical energy in that
convention; it is not a harmless numerical normalization.

The table stores `log(eps + energy_shift)` to permit logarithmic interpolation
when physical eps is negative. Returned energy is
`exp(interpolated_log_energy) - energy_shift`. The interpolation shift must
not be added to evolved tau, stress-energy, enthalpy, or the atmosphere energy.
A genuine change in the physical rest-mass/energy convention would require
consistent changes beyond this interpolation operation.

The removed `eps_min >= 1e-27` floor conflicted with that representation.
Both global inverse endpoints are now bounded before taking the logarithm;
the fixed-rho/Ye temperature-edge bounds are still applied by inversion.
Ideal gas has no interpolation shift and retains its nonnegative energy range.
Startup rejects a supported EOS domain with eps_min <= -1, since the current
GRMHD and recovery bounds require positive `1 + eps` throughout that domain.

## Conservative admissibility and recovery

### Tau

Implementation: [Con2PrimFactory/src/c2p.hxx](Con2PrimFactory/src/c2p.hxx).
Let D and tau be densitized, Btilde the densitized magnetic vector, and
sqrtg the spatial volume factor. The trigger is

```
tau_mag = 0.5 * g_ij Btilde^i Btilde^j / sqrtg
tau < tau_mag + D * min(0, global_physical_eps_min)
```

Only when this fails is tau replaced by

```
tau_mag + D * eps_local_min(rho_est, Ye_est) + sqrtg * tauFluid_atmo
rho_est = clamp(D/sqrtg, EOS rho range)
Ye_est  = clamp(DYe/D, EOS Ye range)
```

This separates a global trigger from the local target at an estimated density,
following the metric overload of FIL's limiter. A local minimum at D/sqrtg
is not itself a bound at every possible recovered rho. The limiter is a
pre-recovery repair, not a sufficient existence theorem for arbitrary GRMHD
states. The existing momentum/total-energy cap remains. It is called only
inside `if (call_c2p)`, so a completed atmosphere reset is not modified again
by that outer pre-C2P limiter.

### Palenzuela

Both trial evaluation and final primitive recovery use the local EOS energy
range. The previous use of `atmo.eps_atmo` as a generic recovered-energy floor
is removed. Ye is assigned before the local energy query. The raw recovered
energy is retained long enough to diagnose non-finite/out-of-domain results.
The old argument name `reject_nonpositive_eps` is retained for compatibility,
but its policy now rejects out-of-local-range energy, not valid negative
physical energy. Final closure/floors use the same shared path as Noble and
RePrimAnd, followed by conservative recomputation when required.

### General-EOS Noble

Reference inspected:
`../grmhd_con2prim/NR_2D_Noble.c`, `eos_quantities_general`.
With x=v^2, W=(1-x)^(-1/2), rho=D/W, and Z=rho*h*W^2, the EOS is asked for
the state at h=Z*(1-x)/rho. Define chi=(dP/drho)_eps and
kappa_e=(dP/deps)_rho (this derivative is unrelated to evolved entropy kappa).
The implementation uses

```
dP/dZ = (kappa_e/rho) * (1-x) / (1 + kappa_e/rho)
dP/dx = [-0.5*D*W*chi - 0.5*kappa_e*(Z + P*W^2)/rho]
        / (1 + kappa_e/rho)
```

The enthalpy inverse is bracketed in the local energy interval, uses Newton
steps when safe and midpoint steps otherwise, and reports convergence.
Noble backtracks trial steps that leave the EOS domain and checks the final
momentum and energy residuals; a tiny shortened step alone is not accepted.

The CompOSE reader does not supply usable values in every pressure-derivative
column. Noble therefore differentiates the actual multilinear interpolants
of log(P) and log(eps+shift), rather than relying on those columns. For table
coordinates R=log(rho), L=log(T), these are

```
dP/deps = P/exp(logeps_shifted) * (dlogP/dL)/(dlogeps_shifted/dL)
dP/drho = P/rho * [dlogP/dR - (dlogP/dL)*(dlogeps_shifted/dR)
                             /(dlogeps_shifted/dL)]
```

Nonpositive/nonfinite energy-temperature derivatives reject this Noble
evaluation, allowing the configured backup policy to handle it. Derivatives
are stencil derivatives, one-sided at table knots; they are not smooth
derivatives of an underlying fitted thermodynamic potential.

### RePrimAnd and shared finalization

RePrimAnd no longer forces the lower enthalpy bound to one. With nonnegative
pressure it uses `1 + min(0, global_eps_min)`, allowing physical states with
h<1. Finalization closes the recovered independent variables, applies the
selected thermal floor and velocity/magnetic policy, and rejects non-finite
results. Temperature-primary magnetic heating is bracketed above the current
T; the caller rejects an unattainable pressure target. The bounded search
does not implement a general multi-root P-to-T inverse, and can reject a
target if the Tmax endpoint is below it even when a nonmonotone interior
pressure maximum exists.

Limited/excised states now reclose thermodynamics and refresh the electric
field after velocity changes. The existing black-hole and magnetic limiter
policies were not redesigned. In extreme states EOS-domain saturation and
magnetization policy can still compete; those cases need dedicated tests.

The entropy C2P class no longer reads Gamma either, but an EOS must implement
its kappa inversion APIs to support that solver. **Tabulated entropy recovery
remains disabled** by validation because those inversions remain unsupported.
EOS-agnostic algorithms do not mean every EOS implements every capability.

## Initial data, face reconstruction, and seeding

- Tabulated ID closes from bounded rho/T/Ye, removes the independent global
  pressure/energy patches, and stores kappa. Atmosphere points use the complete
  canonical thermal state and zero velocity. The subsequent initial
  primitive-to-conservative conversion constructs the conservative fields.
- Ideal-gas primitive checking retains temperature-, pressure-, or
  energy-authority choices. Non-finite or atmosphere-classified primitives
  receive a full reset. The cold initial-data path is retained.
- Saved C2P seeds now derive temperature from saved rho/eps/Ye. The invalidated
  HydroBaseX temperature field is not read to construct the seed.
- Face reconstruction retries the configured lower-order method for invalid
  independent states, then closes or fully resets each face. Zero Ye/T and
  negative physical eps are not rejected solely by sign when allowed by the
  EOS. Evolved kappa and sound speed are derived from the finalized face state.
- HLLE normalizes its weights before multiplication and handles degenerate
  zero/subnormal speed separation without dividing by zero. This is not the
  same dissipation policy as FIL's configurable weak-speed fallback.
- The AsterSeeds beta-floor path stores complete temperature-primary
  thermodynamics and kappa; both scheduled call sites declare their accesses.
  Its existing density-from-pressure/temperature operation and optional
  coorbiting-velocity feature are not a new general-EOS seeding implementation.

When `interpolate_failed_c2p=yes`, failed-cell neighbour inputs and flags are
snapshotted before any repair. Only finite successful neighbours contribute;
conditional sums avoid `0 * NaN`. Averaged rho/eps/Ye are reclosed, local floors
and velocity limits applied, and conservatives including DYe/DEnt rebuilt.
Saved primitives and auxiliary velocity/diagnostic fields are updated as well.
No-neighbour points retain their existing failure handling. Magnetic fields
are preserved; the staggered field is not overwritten. The seven snapshot
grid functions are allocated only when this optional path is enabled. MPI,
AMR, and checkpoint behavior of this path still require execution tests.

## FIL comparison: inspected behavior, not numerical equivalence

This is a source-grounded comparison of the local checkout, not a claim that
an earlier conversation table has been reproduced verbatim or that the two
executables have been benchmarked.

| Area | Local FIL reference | AsterX repair / remaining distinction |
|---|---|---|
| Shifted tabulated energy | `Margherita-EOS/src/3D_Table/tabulated_implementation.hh`: `find_logtemp_from_eps`, energy return routines | Keep physical eps unshifted outside log interpolation |
| Local energy bounds | Same file: `eps_range__rho_ye` evaluates temperature edges at rho/Ye | Use `range_eps_from_rho_ye`, not atmosphere energy |
| Pressure inversion | Same file explicitly rejects generic eps from rho/P | Table pressure-primary reconstruction/floors rejected; ideal gas analytic path retained |
| Conservative tau | `Margherita-EOS/src/c2p/margherita_c2p_mhd_t.hh`: metric overload of `limit_tau_and_return_stilde_sq_max` | Port the global trigger/local-target distinction into AsterX densitized variables; retain AsterX margin |
| Other tau overload | Same file has an array-based overload with a different density estimate/local trigger | Not silently treated as identical; this patch follows the metric overload |
| Atmosphere reset | `IllinoisGRMHD/src/driver_conserv_to_prims.C`: rho/eps/Ye reset and EOS pressure with static physical velocity | Complete shared AsterX reset, preserving AsterX radial grading and Valencia velocity convention |
| Sound speed and flux | `IllinoisGRMHD/src/mhdflux.C` receives left/right csnd2 and constructs HLL speeds | Derive face cs2 once from the same finalized state |
| Weak-speed handling | Same flux routine uses `speed_eps`, replacing sufficiently small speeds by unit bounds | AsterX only gains its explicit degenerate HLLE guard; no wholesale dissipation switch |
| Reconstruction | FIL has `reconstruct_set_of_prims_WENO5.C`; AsterX retains selectable ReconX methods | Fix EOS state consistency; do not claim that changing to FIL's reconstruction alone fixes thermal evolution |
| Diagnostics | FIL C2P tracks error bits and repair statistics | Add optional AsterX event sampling and EOS API/root counters; definitions differ |
| Noble general-EOS Jacobian | Separate `grmhd_con2prim/NR_2D_Noble.c`, not asserted to be FIL's active solver | Use its general derivative relations, with EOSX interpolation/enthalpy closure |

Neither this comparison nor the unit-test definitions establish an explanation
for every observed low-temperature ring. EOS consistency is one necessary
piece; surface advection, truncation error, recovery, resolution, and floor
policies still require controlled evolution comparisons.

## Diagnostics

Both diagnostic switches default to zero. Existing production parfiles are
unchanged. A diagnostic trial can set, for example:

```text
EOSX::eos_call_diagnostics_every = 64
AsterX::repair_every = 64
```

`EOSX::eos_call_diagnostics_every` is a startup setting (`STEERABLE=never`).
When enabled, managed counters accumulate on each MPI rank after EOS setup.
Reporting at POSTSTEP synchronizes device work and prints API categories:
pressure, energy, temperature, physical entropy, evolved kappa, sound speed,
pressure derivatives, local energy-range queries, table inversions, table
root-function evaluations, enthalpy inversions, and enthalpy iterations.
Nested calls count separately: P(T) calling eps(T) and P(eps) is not one unique
thermodynamic state evaluation. Chemical-potential/composition queries and
all possible EOS entry points are not exhaustively instrumented. Table loading
is excluded; enabled startup tests may contribute to subsequent counts.
No global MPI reduction or checkpoint persistence of counters is provided.

`AsterX::repair_every` samples C2P/flux invocations on selected iterations.
Each invocation owns its counter storage. Output includes rank, iteration,
routine label, and flux direction; multiple RK stages/levels/tiles can produce
multiple records. C2P loops include ghost points, so these are **work/event
counts, not unique physical-cell counts or mass budgets**.

| Output | Meaning / coverage |
|---|---|
| `cell_atmo` | Atmosphere decisions in the sampled main C2P path |
| `face_atmo` | Individual left/right thermal face resets |
| `rho_clamp`, `T_clamp`, `Ye_clamp`, `eps_clamp` | Instrumented final recovery/face closures, including attempted solver closure; not every internal table endpoint evaluation |
| `tau_repair` | Main pre-C2P limiter changed tau |
| `primary_fail` | Primary recovery failure, including a deliberately disabled primary |
| `backup_call`, `backup_fail` | Configured second-solver calls/failures; not a separate entropy-backup summary |
| `cons_recompute` | At most one recorded rebuilt/reset decision per main C2P point; not a count of every `from_prim` invocation |
| `loworder_face` | Interface switched to its configured lower-order reconstruction |
| `pplim` | Sampled limiter blend had theta<1; limiter algorithm unchanged |

Initialization, beta-floor, optional neighbour repair, every magnetic/BH
limiter event, and every pressure-face energy clamp are not separately counted
by the main repair sampler. EOS calls inside these paths still contribute to
enabled EOS API counters where instrumented. These limitations matter when
interpreting totals. Atomics, allocation/copying on sampled calls, and report
synchronization have overhead; measure performance with diagnostics off too.

## Tests added inside the relevant thorns

| Thorn / source | Coverage written |
|---|---|
| `EOSX/src/test.cxx` | Uniform/nonuniform multilinear value and derivative checks; ideal-gas closure for three Gammas including zero energy; bounds, enthalpy failure, API-counter smoke checks; 64 active-EOS device round trips |
| `Con2PrimFactory/src/test_eos_repairs.cxx` | Synthetic shifted table with zero derivative columns and negative physical eps; EOS/enthalpy round trips; analytic and finite-difference Jacobians; independent Noble/Palenzuela/RePrimAnd recovery inputs across density, temperature, velocity and magnetic samples; atmosphere fixed point; cold/thermal grading; non-unit-metric tau tests; bounded pressure-floor heating |
| `ReconX/src/test.cxx` | Constant (including negative), linear, reflected and jump profiles for minmod, MC, WENOZ, WENOZp, MP5 and PPM; finite-output and TVD checks |
| `AsterX/src/unit_tests/hlle.cxx` | Zero-speed, one-sided upwind and symmetric-speed flux limits |
| `AsterX/src/unit_tests/thermo.cxx` | Distinct left/right cold atmosphere grading, pressure face closure, acoustic speed consistency |
| `AsterSeeds/src/test.cxx` | Closure and evolved kappa after density changes at fixed T/Ye; ideal gas and active table |

The synthetic table is `eps=T-0.05`, `P=rho*T`, with a thermodynamically
consistent simple sound speed. It deliberately has negative physical energy
and h<1 at low T. Each energy solver receives an independent conservative
input, preventing one solver's repair from making the next solver's test pass.
This synthetic EOS is a mathematical regression fixture, not a nuclear EOS.

The standalone startup input is
[AsterX/test/eos_repairs.par](AsterX/test/eos_repairs.par). It stops at iteration
zero, activates all five suites, and needs no external table for the synthetic
tests. It has not been run and has no generated reference-output directory.

For actual DD2 coverage, add these switches to a **copy** of an otherwise valid
tabulated run input after rebuilding:

```text
EOSX::unit_test = yes
Con2PrimFactory::unit_test = yes
ReconX::unit_test = yes
AsterX::unit_test = yes
AsterSeeds::unit_test = yes
```

The active-EOS tests then query the loaded table. Do not replace the ideal-gas
test input's EOS name alone: its Balsara initial-data routine is ideal-gas-only.
The real table needs its normal table filename/format, initial data, and mode
settings. Tests abort on failure rather than merely printing discrepancies.

Not yet covered by executed tests: everything above. Additional integration
coverage still required includes real DD2 C2P parameter sweeps; both HDF5
reader formats; table-axis and monotonicity assumptions; low-T inversion
conditioning; extreme magnetization/velocity; entropy fallback for ideal gas;
full ePPM/Godunov grid-function paths; scheduled seeding; forced failed-C2P
neighbours; MPI/AMR/subcycling; checkpoint/restart; and stellar evolution.
The production table class still uses the existing uniform-axis interpolator;
testing the generic nonuniform interpolator does not change that assumption.

## Static review and handoff

Reviewed changed source/schedule paths, energy conventions, local/global
limits, dispatch, declarations, and new test registrations. Follow-up commits
record review findings rather than hiding them in a squashed patch: frozen
neighbours, missing namespace import, cold-EOS dispatch, invalidated seed
temperature, stale electric field after limiting, parameter steering guards,
and counter interpretation.

`git diff --check origin/dev` passes. Source search finds no
`GammaIdealFluid`, direct `eos_3p->gamma`, or `eos_3p->gm1` dependency in the
C2P files, and no shortened `eos_call_every` parameter. No compiler, linker,
Cactus schedule-generation, or runtime validation has been performed.

Suggested validation order after user review/rebuild:

1. Run the startup test input, then the active-DD2 startup suites.
2. Run static atmosphere tests with zero and nonzero grading, both EOSs,
   unmagnetized and magnetized, at cell and face locations.
3. Exercise each energy C2P independently on the same real-table states;
   compare primitive/conservative residuals and count fallback frequency.
4. Run short controlled star tests with the original face cutoff retained,
   monitoring temperature, baryon mass, constraints, resets, tau repairs,
   and API/root counts. Compare like-for-like resolution and output times.
5. Repeat MPI/AMR/subcycling and restart cases; force the optional neighbour
   path and scheduled seeding separately.
6. Only then compare `recon_use_atmo_tol=yes` and PPLIM/reconstruction choices
   as isolated experiments. Do not combine them with the first repair test.

These repairs intentionally change formerly inconsistent floor/surface states
and the tabulated Noble calculation. They cannot promise unchanged evolution
or an improved temperature profile before those comparisons. No movies or
plots were regenerated, and no simulation was submitted for this series.
