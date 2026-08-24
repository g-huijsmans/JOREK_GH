# Historical JOREK-GK Poisson validation

August 2026

## Purpose and references

This validation compared model011's direct physical-potential Galerkin
Poisson implementation with the historical JOREK-GK particle
implementation. It covered the isolated weak form, particle initialization,
deposition, the full particle application, and cross-solves using identical
final right-hand sides (RHSs).

The historical source state was
`/pitagora_work/FUPA2_MHDnoELM/ghuijsma/marconi_iter_irene/jorek`. That tree
already contained significant local modifications, and those modifications
were part of the executable being validated. They were not created by this
validation. Consequently, the comparison used an isolated "Route B" copy of
the actual working-tree source state rather than a clean historical commit.
The temporary rebuild was byte-for-byte identical to the historical
executable, whose SHA256 was:

```text
576b8713ce5a2deb3a25564df5b5bce55292fc92d0d4de75812359a94b324629
```

Related permanent references are:

- [the original coefficient comparison](poisson_original_comparison.md);
- [the harmonic and cosine/sine regression](test_poisson_harmonics.f90);
- [the solver-cycle, family, and direct-RHS regression](test_poisson_solver_cycle.f90);
- [the restored particle projection diagnostics](../../../particles/tests/test_projection_diagnostics.f90).

## Isolated matrix and weak-form validation

The first comparison used coefficient-matched inputs and excluded the full
particle application. After correcting discrepancies in the model011 kernel,
the historical and new matrices agreed approximately at roundoff:

| Comparison | Relative error |
|---|---:|
| Filters disabled | `~1.8e-15` |
| Individual filter terms | `<~3.4e-15` |
| Full filtered matrix | `~2.66e-15` |
| Potential | `~1e-15` |
| Controlled direct RHS | `0` |

The substantive corrections were:

- full-field polarization subtraction;
- same-harmonic cosine/sine coupling;
- correct coefficients for the `n=0` and nonzero-harmonic filters;
- physical Fourier normalization of `2*pi` for `n=0` and `pi` for `n>0`.

This established agreement for the isolated coefficient-matched weak form;
it did not yet validate the complete particle application.

## Particle initialization and robustness

The full application comparison exposed two particle-temperature differences.
First, `T_particles` support had disappeared from
`initialise_particles_H_mu_psi`. It was restored with the historical
precedence:

1. use `T_maxwell` when supplied;
2. otherwise use `T_particles(psi)` when supplied;
3. otherwise interpolate the physical `var_T` field.

Second, `T_ions` and `T_electrons` initially used the normalization for the
single-temperature JOREK field. JOREK uses

```text
T = Ti + Te
```

for that field, while `T_ions` and `T_electrons` return individual species
temperatures. The callbacks therefore use the individual-species conversion

```text
1 / (EL_CHG * MU_ZERO * central_density * 1e20)
```

without the factor `1/2` present in the single-temperature normalization.
After this correction, historical and current thermal moments agreed closely.

The historical run reserved no inactive markers (`zero_fraction=0`), while
CURRENT production normally reserves ten percent (`zero_fraction=0.1`). For
the controlled comparison only, CURRENT used zero. This produced identical
active-particle counts and essentially particle-for-particle agreement. A
representative matched-profile case contained 987242 active ions and 987497
active electrons; position and weight distributions agreed near machine
precision. Production retains `zero_fraction=0.1`.

Two independent robustness bugs were also corrected:

- invalid gyro points could leave negative element indices entering
  interpolation; invalid points are now skipped during gyro averaging;
- the guiding-centre RK4 lost-particle element state was not propagated from
  intermediate searches; the relevant RK stages now propagate `i_elm`.

These fixes were required for robust million-particle runs, but they were not
the ultimate cause of the potential discrepancy. The validation also restored
`proj_vpar`, `proj_Pressure`, and `set_vtk_active` behavior.

## Deposition and direct-RHS validation

After matching the particle populations, the apparent historical `Te_eV/zn0`
factor was shown not to be a missing normalization in CURRENT. The historical
formulation applies the same local `T_e/n` transformation to its RHS and its
polarization operator.

The first common-convention deposition difference was the known time
normalization

```text
t_old / t_current
  = sqrt(MASS_PROTON / ATOMIC_MASS_UNIT)
  = 1.003631603838641.
```

After applying this analytic factor, without fitting vectors, the deposition
differences were approximately:

| Contribution | Relative difference |
|---|---:|
| Ion | `~6e-9` |
| Electron | `~1e-9` |
| Net ion plus electron | `~5e-8` |

The CURRENT direct RHS also agreed with a full-mesh projection-reference path
to approximately `4.7e-9` relative on 135245 production degrees of freedom.
This validated CURRENT particle deposition, global scatter, and direct Poisson
RHS construction.

## Constant-equilibrium-`T_e/n` experiment

An artificial equilibrium was constructed with complete Hermite fields

```text
T_h = C rho_h
```

for every value and derivative degree of freedom. Thus `T_e/n` was constant to
roundoff. Particle phase space was held independently at
`T_maxwell = 500 eV`, so the artificial equilibrium temperature changed the
weak form but not particle initialization.

The potential correlation changed from approximately `0.76-0.81` with the
real variable profiles (depending on the exact matched setup) to `0.9999164`
with constant equilibrium `T_e/n`. The latter full-field comparison gave:

| Metric | Value |
|---|---:|
| Relative L2 | `7.76e-3` |
| Maximum difference | `2.17e-4` |
| Correlation | `0.9999164` |

This localized the large historical discrepancy to the spatial `T_e/n`
transformed-test formulation.

## Historical 10 eV floor defect

The historical implementation evaluated approximately

```text
T0_g = max(T0_g * Tev_norm, 10.d0)
```

but retained derivatives of the raw field:

```text
T0_s = T0_s * Tev_norm
T0_t = T0_t * Tev_norm.
```

Where `T_raw < 10 eV`, it therefore used `T_eff = 10 eV` together with
`grad(T_eff) = grad(T_raw)` in

```text
grad ln(T_e/n) = grad(T_e)/T_e - grad(n)/n.
```

That derivative is inconsistent with the clipped temperature. In the matched
case, 97935 of 535600 Gaussian points (18.29 percent) were clipped. The minimum
raw temperature was 1.5521 eV, and no raw temperatures were negative.

Lowering the floor enough to make it inactive made agreement much worse, so
the floor itself should not simply be removed. The surgical test retained

```text
T_eff = max(T_raw, 10 eV)
```

but set `grad(T_eff)=0` where the floor was active. It changed the
OLD-versus-CURRENT comparison to approximately:

| Metric | Corrected-floor value |
|---|---:|
| Relative L2 | `0.081547` |
| Maximum difference | `4.55e-4` |
| Correlation | `0.998935` |

The corresponding historical correlation was about `0.755` in the matched
constant-particle-temperature, variable-profile case. The dominant historical
implementation defect was therefore the inconsistent derivative of the
clipped temperature, not the mere existence of the 10 eV floor.

## Corrected Petrov formulation versus CURRENT Galerkin

The corrected historical transformed-test formulation is

```text
A_P(v, phi) = integral R * (m_i T_e)/(e B^2)
              * [grad(v).grad(phi)
                 + v grad ln(T_e/n).grad(phi)] dR dZ.
```

The production CURRENT formulation is

```text
A_G(v, phi) = integral R * (m_i n)/(e B^2)
              * grad(v).grad(phi) dR dZ.
```

Historical testing effectively multiplies the test function by
`f=T_e/n`. The two forms can represent the same strong equation at the
continuous level. In a finite-dimensional space, however, `f v_h` is
generally not in the original test space. The historical method is therefore
a transformed-test, or Petrov-Galerkin, discretization, while CURRENT directly
discretizes the physical-potential Galerkin weak form. This is a legitimate
discretization distinction, not evidence that the historical continuous
formulation is mathematically wrong.

## CURRENT_PETROV cross-implementation test

For diagnostic purposes, the corrected historical Petrov form was
reimplemented using CURRENT basis functions, quadrature, assembly, boundary
conditions, and direct solver. Across all 33475 elements, corrected OLD and
CURRENT_PETROV element matrices agreed with relative Frobenius error
`1.80239e-8`. Agreement was uniform across value-value, value-derivative,
derivative-value, and derivative-derivative blocks. The remaining element
difference was traced to historical deuterium-mass precision, not to a
different integrand.

This establishes that CURRENT infrastructure can reproduce the corrected
historical Petrov operator.

## Common-RHS cross-solves

Independent corrected-OLD and CURRENT_PETROV solves initially differed by
relative potential L2 `8.70845e-3`, with correlation `0.9998835`. Their small
net RHSs result from strong ion/electron cancellation and are consequently
very sensitive to tiny independent particle-state differences.

Both implementations were then supplied with exactly the same final
coefficient vector immediately before their direct solves:

| Canonical RHS | Relative potential L2 | Correlation |
|---|---:|---:|
| CURRENT RHS in both systems | `3.8905e-5` | `0.9999999977` |
| OLD RHS in both systems | `2.8994e-5` | `0.9999999984` |

More than 99.5 percent of the apparent 0.87 percent independent-solve
difference therefore came from independently accumulated, strongly
cancelling particle RHS vectors. With a common RHS, corrected OLD Petrov and
CURRENT_PETROV reproduce one another to the expected numerical accuracy.

## Conclusions

1. The dominant original full-application discrepancy was caused by the
   historical inconsistency between clipped temperature values and unclipped
   derivatives.
2. Particle initialization and deposition differences were independently
   identified and controlled.
3. CURRENT direct deposition, global scatter, and Poisson RHS construction
   were validated.
4. Corrected OLD Petrov and CURRENT_PETROV agree tightly when supplied the
   same RHS.
5. The remaining CURRENT_PETROV-versus-CURRENT_GALERKIN relative L2 difference
   of approximately 8.19 percent is the genuine finite-dimensional
   transformed-test Petrov-Galerkin versus direct Galerkin difference.
6. There is no remaining evidence of a model011 Poisson-operator
   implementation bug.
7. CURRENT retains the direct physical-potential Galerkin formulation for
   production.

## Limitations

- The full-application validation discussed here primarily covers
  `n_tor=1`, or `n=0`.
- Nonzero-harmonic cosine/sine coupling is covered separately by the linked
  harmonic and solver-cycle regression tests.
- Direct mode-family structure currently requires contiguous component
  intervals.
- The proton-mass versus atomic-mass-unit time-normalization convention remains
  a small uniform difference unless explicitly aligned.
- Strongly cancelling particle RHS vectors are highly sensitive to tiny
  independent particle-state differences.
