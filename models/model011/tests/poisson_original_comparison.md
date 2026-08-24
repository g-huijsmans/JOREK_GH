# Original gyrokinetic Poisson comparison

## Reference and test scope

The authoritative source was inspected read-only at
`/pitagora_work/FUPA2_MHDnoELM/ghuijsma/marconi_iter_irene/jorek`.
The isolated snapshot and matched input copies are in
`/tmp/jorek_poisson_reference_test`: the copied inputs used for configuration
comparison are `reference_input` and `current_input`, and the executable
numeric oracle is `test_poisson_original_reference`.

`particles/examples/testing.f90` does not use the fluid time-stepper for this
solve. It constructs a `projection` with ion polarization enabled and calls
`with(sim,jorek_feedback)`. Consequently the active reference implementation
is `particles/diagnostics/mod_project_particles.f90`:

- `prepare_mumps_par`: nonzero-mode matrix and boundary rows;
- `prepare_mumps_par_n0`: n=0 matrix and boundary rows;
- `project_only`: element-load scatter, MPI reduction, MUMPS solve, and field storage;
- `loop_particle_gc_local` in `testing.f90`: particle feedback deposition.

The copied full program was not built or run. It pushes particles before the
solve and depends on a large machine-specific application state, so it is not
a clean strict Poisson oracle. The validation instead uses the allowed common
saved-FE-RHS fallback and an independent plane-sampled implementation of the
active integrand in `test_poisson_harmonics.f90`.

## Inputs and disabled physics

The supplied old input has one scalar `n_particles`; `testing.f90` creates the
ion and electron groups, assigns charges +1/-1, assigns deuterium mass and a
100-times-lighter electron mass, initializes both with Sobol sequences, and
normalizes statistical weights with `adjust_particle_weights`. The supplied
new input declares those fields per `part_group_configs(i)` and adds
`nsubstep_electrons`.

The copied inputs match one step, `tstep_particles=2e-8`, 1,000,000 particles
per species, ion mass 2.01410178 u, electron mass 0.0201410178 u,
`central_density=0.8`, `central_mass=2.01410178`, and the old example's filter
values. The old program itself fixes all four binary collision switches and
both small-angle collision switches to false. It also fixes
`fuelling_rate=0`, `heating_power=0`, and `V_mom_ion_source=0`. Those are the
actual names and values; no guessed namelist switches were added.

| Physical item | old/reference declaration | current/new declaration |
|---|---|---|
| species | two groups created by `testing.f90` | `part_group_configs(1:2)%id` |
| particle type | allocated as `particle_gc_vpar` in code | `%type="particle_gc_vpar"` |
| charge | +1 ion, -1 electron passed by the two initializer calls | `%Z=1`; the group role/initializer supplies the particle sign |
| ion mass | `atomic_weights(-2)` in code | `%mass=2.01410178` |
| electron mass | `atomic_weights(-2)/mass_ratio`, `mass_ratio=100` | `%mass=0.0201410178` |
| particle count | scalar `n_particles`, split into both groups | per-group `%n_particles` |
| density/weight | density-profile integral followed by `adjust_particle_weights` | same physical profile/weight adjustment path |
| temperature | `T_particles=T_ions/T_electrons` callbacks in code | initializer callbacks; no per-species temperature in the supplied block |
| initialization | hard-coded `sobseq_rng`, uniform-space rejection | current initializer path plus explicit group metadata |
| substeps | `nsubstep_particles`; electron setting absent/commented | adds `nsubstep_electrons` |
| sources/collisions | hard-coded source rates and logical collision parameters | no corresponding enabled option in the copied controlled input |

Toroidal layout is compile-time configuration, not an input field. The
authoritative reference snapshot currently has `n_tor=1`, `n_period=1`, and
`n_plane=1`, while the current tree has 3, 1, and 4. Those application builds
therefore cannot be mode-matched by editing the copied namelists alone. The
isolated oracle explicitly uses `n_period=1`, 64 sampling planes, and modes 0,
2, and 6; it exercises the reference nonzero-mode routine without modifying
the reference compile-time settings.

The original Sobol initialization is deterministic within that executable,
but it is not guaranteed to produce the same particles in the current code.
No comparison below uses independently initialized particle sets.

## Active weak form

Let `L f = f_RR + f_R/R + f_ZZ`,
`D f = (F0 f_phi/R + f_R psi_Z - f_Z psi_R)/R`, and
`B2=(F0^2+psi_R^2+psi_Z^2)/R^2`. For a test function `v` and trial function
`p`, the active nonzero-mode integrand, before `R J w`, is

```
filter_perp grad(v).grad(p)
+ filter_hyper L(v)L(p)
+ filter_par D(v)D(p)/B2
+ Fpol [(grad(v)+v grad(log(T/n))).grad(p) + v_phi p_phi/R^2]
- Fpol D(v)D(p)/B2
```

where `Fpol=m_i T_eV/(e B2)`. Adiabatic-electron and Pade terms are disabled
by `testing.f90`. For n=0, the polarization parallel subtraction is commented
out, there is no toroidal derivative, and an active hard-coded contribution
adds `D(v)D(p)/B2` where normalized psi is below 0.64. The configured n=0
perpendicular, hyper, and parallel filters are otherwise analogous.

The direct reference RHS is scattered from element-local particle loads
without a projection mass solve. Each particle contribution contains its
charge, statistical weight, basis value, real Fourier value, `t_norm/F0`, and
the reference variable's local `T_eV/n`. `project_only` sums shared nodes and
MPI ranks and passes that vector directly to MUMPS. The current physical-load
representation instead deposits charge times `t_norm/F0` and divides the
assembled vector by `central_density*1e20`. A strict comparison therefore uses
one common assembled RHS in the current physical normalization; the two raw
particle-deposition arrays are not the same algebraic variable.

The reference matrix multiplies the plane sum by `2*pi/n_plane`, giving
physical Fourier norms `2*pi` for n=0 and `pi` for each nonzero sine/cosine
component. Direct particle delta loads have no compensating plane sum.

## Analytic harmonic form

For a nonzero physical mode `k`, define

```
d_i = (v_i,R psi_Z - v_i,Z psi_R)/R
a   = F0/R^2
q   = (filter_par - alpha)/B2
```

where `alpha` is the polarization coefficient in the current physical
normalization. The corrected real harmonic blocks are

```
A_ij = RJ [(alpha+filter_perp) grad(v_i).grad(v_j)
           + filter_hyper L(v_i)L(v_j) + q d_i d_j]
B_ij = RJ [alpha v_i v_j/R^2 + q a^2 v_i v_j]
C_ij = RJ q a (d_i v_j - v_i d_j)

K_cc = K_ss = pi (A + k^2 B)
K_cs =  pi k C
K_sc = -pi k C
```

`C` is skew-symmetric, so the complete real matrix is symmetric. Axisymmetry
does not make this cosine-sine block zero: the toroidal component of
`B.grad` rotates sine into cosine. Different physical mode numbers remain
uncoupled. For n=0, the block is `2*pi A0`, uses the n=0 filters, includes the
active core parallel coefficient, and has no polarization subtraction.

The original `T/n` variable and the current density-normalized physical
potential are not coefficient-identical. A constant scalar filter conversion
exists only when their normalization ratio is constant: after scaling the
reference equation by `c=T_eV/rho`, `filter_current=filter_original/c`.
With spatially varying `c`, no single namelist conversion reproduces the
reference filter strength everywhere. The numerical kernel comparison below
matches the polarization coefficient explicitly and tests the geometric and
Fourier operator independently of that variable-definition issue.

## Boundary conditions

The reference nonzero-mode solve applies homogeneous Dirichlet rows to
boundary types 1, 2, 3, 4, 5, and 9. With the flags used by `testing.f90`, n=0
applies them to types 1, 2, and 3, but not 4, 5, or 9. The current generic
boundary code follows `bcs(boundary)%dirichlet%u` and can apply all of them,
and also has `keep_n0_const` and explicit axis-basis transformation. Thus a
full-mesh run with `bcs(:)%dirichlet%u=.true.` does not exactly match the
reference n=0 boundary selection. The controlled matrix/solution test applies
the same homogeneous Dirichlet DOFs to both matrices.

## Numerical results

The independent reference sums 64 toroidal planes and three poloidal
quadrature samples for modes 0, 2 cosine/sine, and 6 cosine/sine. Reported
columns are maximum absolute and global maximum-relative errors:

| case | matrix abs | matrix rel |
|---|---:|---:|
| filters off | 4.62e-14 | 1.79e-15 |
| perpendicular only | 2.31e-14 | 8.69e-16 |
| hyper only | 2.49e-14 | 9.63e-16 |
| parallel only | 3.20e-14 | 9.41e-16 |
| active n=0 core term only | 4.62e-14 | 1.79e-15 |
| all terms | 9.24e-14 | 2.66e-15 |

For a common deterministic RHS and identical Dirichlet elimination, the dense
comparison solve gives max absolute 1.53e-16, max relative 1.09e-15, and L2
2.0e-16. The physical projected-RHS quadrature gives max absolute 4.81e-34
and relative 1.09e-15; the saved direct FE RHS normalization is exact for n=0,
cosine, and sine. This small element is also used as the mapped global system;
a full application sparse-matrix dump was not produced.

The two-rank solver-cycle regression passes and retains factorization,
`solve_only`, direct-vs-projected RHS, and absolute-phi storage checks. Its
changed-RHS response is 3.06e-13; repeated solve, linear scaling, direct RHS,
storage, and untouched-field errors are zero.

## Corrections and performance

The comparison demonstrated and corrected four active discrepancies in the
current analytic path:

1. missing full-field parallel subtraction from ion polarization;
2. missing sine-cosine coupling generated by that field-aligned derivative;
3. use of n=0 filter coefficients for every harmonic and omission of the
   active n=0 core parallel term;
4. plane-sum (`n_plane`, `n_plane/2`) instead of physical (`2*pi`, `pi`)
   Fourier normalization for the direct-load equation.

The analytic construction evaluates equilibrium and basis data once per
poloidal quadrature point and scatters four reusable poloidal blocks. The
reference evaluates the full integrand on every toroidal plane. Its dominant
element work scales with `n_plane` times the local basis-pair work (plus
Fourier sampling), whereas the analytic work is one block build plus linear
scatter per retained harmonic. No wall-clock speedup is quoted because the
full copied reference executable was intentionally not built.
