# Coverage of the supplied draft

The numbering below follows `inputs/draft.pdf`: Table 3 contains the
accessory families, Table 4 the selected CM images, Table 5 the assumed
initial coefficients, Table 6 the expansions, and Table 7 the certified
coordinates. Some names inside the supplied Nelson source retain earlier
table numbers; the package's coverage uses the supplied PDF.

All executed checks pass. Exact calculations use integer, rational,
algebraic-number, or formal power-series arithmetic. Numerical comparisons
are labeled separately. Statements imported from the literature remain
explicit mathematical inputs.

| Draft claim or computation | Certificate and evidence | Scope of the executed check |
| --- | --- | --- |
| Table 1 models and eigendifferentials; quotient coordinates | `certificates/families.py`; `data/families_certificate.json` | All nine models are squarefree; the stated involutions preserve the models and coordinates and act by signs on the differentials. H is derived from each model, differential, and quotient coordinate. Identification with the Shimura curves uses the cited model theorems. |
| Theorems 1.1-1.2 and equations (1.5)-(1.6), as used for the twelve forms | `certificates/coefficients.py`; `results/coefficients.json` | Implements the exact recursion and checks both original undivided equations and `(u')^2 = L A^2 H(u)` through 30 residual coefficients, using 32 computed terms. The general Schwarzian and descent proofs are the mathematical arguments in the draft. |
| Table 2, Section 4, and Theorem 1.3: recovery of R, L, u_0, u_1 | `certificates/notebooks.py`; both executed notebooks in `results/`; `results/recovered_parameters.json` | Every code cell is freshly executed. The known Schwarzian transformations for 14a, 15a, 21a and parameter recovery for the other nine forms reproduce the table. These computations assume the stated initial coefficients. |
| Table 3: all accessory families and parameter counts | `certificates/families.py`; `data/families_certificate.json`; `docs/families.md` | Derives all twelve families from local Schwarzian principal parts and the descended divisor. Checks local ramification, differential zeros, independent accessory bases, dimensions, and exact equality with the displayed families. Elliptic fiber identifications come from Appendix C and cited geometry. |
| Appendix A, Table 5: 43 stated initial coefficients | `certificates/numerical.py`; `results/numerical.json` | Nelson numerical evaluation at both lattice cutoffs, independent of the recurrence's generation of these values. At cutoff 196 the maximum observed discrepancy is below `1.21e-57`. Exact rational recognition remains a hypothesis. |
| Appendix B, Table 6: all nine tau values | `certificates/nelson_inputs.py`; `results/nelson_inputs.json` | Exact quaternion and maximal-order checks; CM-order discriminants and class numbers; exact fixed-point equations, equality with the displayed radical expressions, and positive imaginary parts. |
| Appendix B, Table 6: twelve Theta values | `certificates/numerical.py`; `results/numerical.json` | At cutoff 196 the maximum observed difference for Theta is below `4.10e-68`; for Theta/(pi Omega_D^2) it is below `6.17e-66`. The eta/Gamma period comparison is also recomputed at high precision. These are numerical agreements. |
| Appendix B: displayed expansions and arbitrary further coefficients | `certificates/coefficients.py`; `results/coefficients_first20.md`; `results/coefficients.json` | All printed coefficients are checked; exactly 20 normalized coefficients per form are displayed. Public functions return any requested count, appendix metadata, or coefficients after multiplying by `(1-w)^2`. Results are exact consequences of the recurrence data. |
| Appendix C, Table 7: all 24 coordinate rows and reciprocal conventions | `certificates/appendix_c.py`; `results/appendix_c.json`; `docs/appendix_c.md` | Explicit published-coordinate changes, model and quotient identities, finite/infinite value substitution, and all genus-two reciprocal pairs. Published CM-value identifications are recorded as theorem inputs. |
| Appendix C, equation (C.1): the Lin-Yang parametrization argument | `certificates/appendix_c.py`; `docs/appendix_c.md` | Full weighted-invariant polynomials are substituted at the three centers; all four coordinates, scaling factors, and equality with the discriminant-6 parametrization are checked exactly. Errthum's two norm values determine the sign. Reduced quadratic forms verify class number one; the moduli/commutant argument is documented with its literature inputs. |
| Consistency of the independently supplied and generated data | `certificates/integration.py`; `results/integration.json` | All twelve notebook-recovered parameter sets match the recurrence data, derived family, Nelson inputs, and Appendix C center. Exact tau and scale expressions are compared across modules. Missing prerequisite records are not replaced with manuscript answers. |
| Section 5, Table 4: four printed polynomial rows | `certificates/cm_claims.py`; `results/cm_claims.json`; `docs/cm_claims.md` | Exact primitive irreducibility, degree/class-number agreement, CM admissibility and selected quaternion input records, plus positive-norm Atkin-Lehner normalizers. The four polynomials remain the manuscript's stated numerical image claims; analytic image identification was not rerun. |
| Section 5: twelve-under-200 source package | `certificates/cm_sanity.sage`; `results/cm_sanity.json` | Small 20-term setup checks for all twelve forms at 128 bits, with selected CM input and reduction checks. The full discriminant scan, analytic integration, recognition, and complete-table fingerprint reproduction were intentionally not run. |

## Numerical values in the recorded run

The run uses 90 decimal digits, an 80-digit weight-function target, and
squared cutoffs 144 and 196. At the final cutoff the maxima are:

```text
initial coefficients:       1.20726833260848333122531058...e-57
Theta:                      4.09806538577807860636703959...e-68
Theta/(pi Omega_D^2):        6.16957959921282684727...e-66
eta/Gamma period comparison: 6.488142e-89
```

The worst coefficient discrepancy when both cutoffs are included is
`1.43265316802199958509...e-47`; the maximum change between cutoffs is
`1.43265316790127275183...e-47`. Thus the saved lower-cutoff run already
satisfies the requested `1e-40` coefficient comparison, while the final
cutoff gives the stronger agreement above. These observations also satisfy
the draft's `1e-45` coefficient and `1e-50` normalized-scale comparisons.

The results quantify numerical agreement with the proposed exact values.
They do not enclose the omitted lattice tail, quadrature remainder, or
rounding through coefficient ratios. See `docs/numerical_scope.md` for the
precise distinction from a rigorous analytic bound.

## Mathematical and execution boundaries

The code certifies the computations used in the stated conditional results;
it does not replace the draft's general proofs or reprove the cited model,
moduli, CM-value, and transitivity theorems. Appendix C includes the
geometric reasoning needed to connect its exact algebra to those inputs.

The Section 5 execution follows the requested sanity-only scope. The
source's full-table regeneration remains available, but its numerical
images and full-table hash are not newly certified by this run. The
approximate 4400-term recognition statement in the draft was not benchmarked
by the 20-term sanity checks.

The original notebooks and external files were preserved. Freshly executed
copies and exported parameter records are in `results/`. The local copy of
the under-200 source has only the documented Sage 9.5 compatibility
adaptations needed by the sanity run; the original archive is unchanged.
