# Computational certificate for the recurrence-formula draft

This package collects the computations used in the supplied draft,
[included here](inputs/draft.pdf). It covers all twelve eigenforms, the
parameter recovery in Section 4, all Table 3 families, the CM expansion
coefficients and scales, and every Appendix C coordinate. Section 5's
twelve-under-200 source is included and has been run only in sanity mode,
as requested.

Start with [the execution report](REPORT.html),
[the claim-to-certificate map](docs/coverage.md), or
[the first 20 coefficients of every form](results/coefficients_first20.md).
The machine-readable overview is [results/summary.json](results/summary.json).
The complete execution log is [results/execution.log](results/execution.log);
[data/source_manifest.json](data/source_manifest.json) records the supplied
file hashes and the small source compatibility adaptations.

## Reproduce the package

From this directory, with SageMath available:

```console
sage -python run_all.py
```

This executes the exact certificates, both notebooks from fresh Sage kernels,
the Nelson numerical comparisons at both cutoffs, the integration checks,
and the small CM sanity checks. It never invokes the full under-200 batch.
The numerical comparison stage took about six minutes in the recorded run;
runtime depends on the machine.

To check the saved report records and their integrity without rerunning the
computations:

```console
sage -python run_all.py --check-recorded
```

The package was executed with SageMath 9.5 on Ubuntu. Use an installed Sage
environment with NumPy, SymPy, mpmath, nbformat, nbclient, and jupyter_client.
The notebook runner creates its own temporary kernel registration bound to
the Sage Python that launches it. No network access is needed for the
computations. Do not use Python's `-O` option: some supplied checks use
assertions.

## Numerical agreement obtained

The recorded run used 90 decimal digits, an 80-digit target for the weight
function, and squared lattice cutoffs 144 and 196. At the final cutoff 196,
the largest observed differences over all twelve forms were:

| Comparison | Observed maximum, rounded upward |
| --- | ---: |
| 43 initial coefficients versus the stated rational values | `1.21e-57` |
| Theta versus the stated period times algebraic value | `4.10e-68` |
| Theta/(pi Omega_D^2) versus the stated algebraic value | `6.17e-66` |

All requested `1e-40` comparisons pass. The initial-coefficient differences
are also below `1e-45`, and the normalized-scale differences are below
`1e-50`. The eta/Gamma period comparison has largest observed difference
`6.488142e-89`. Full values and the changes between cutoffs are recorded in
[results/numerical.json](results/numerical.json).

These are observed numerical differences, rather than rigorous enclosures
of the infinite analytic sums. Exact algebraic checks and numerical evidence
are distinguished in [the numerical scope note](docs/numerical_scope.md).
The recurrence results retain the draft's assumptions on the initial
coefficients and algebraic scale.

## Compute any requested number of coefficients

The normalization is exactly the one in Appendix B:

\[
f_N(z)=(1-w_N)^2\sum_{n\ge0}\frac{c_{N,n}}{n!}(\Theta_Nw_N)^n,
\qquad w_N=\frac{z-\tau_N}{z-\overline{\tau_N}}.
\]

From a Python session launched in the package directory:

```python
from certificates.coefficients import first_coefficients, appendix_expansion, w_coefficients

c = first_coefficients("46a", 20)       # exact Fractions c_0,...,c_19
expansion = appendix_expansion("46a", 20)  # tau, scale, c_n and c_n/n!
w_terms = w_coefficients("46a", 20)     # includes the (1-w)^2 factor
```

Replace `20` by the desired count. Accepted labels are
`14a, 15a, 21a, 26a, 26b, 33a, 34a, 38a, 38b, 46a, 58a, 58b`.
The displayed lists contain 20 coefficients, indexed 0 through 19. The
saved verification computed 32 terms per form and checked both undivided
equations and the first integral through 30 residual coefficients.
See [the coefficient documentation](docs/coefficients.md).

## Individual computations

All commands below are run from the package directory:

| Computation | Command | Main output |
| --- | --- | --- |
| General Schwarzian and recurrence identities | `sage -python certificates/general_identities.py` | `results/general_identities.json` |
| All H functions and Table 3 families | `sage -python certificates/families.py` | `data/families_certificate.json` |
| All Appendix C coordinates | `sage -python certificates/appendix_c.py` | `results/appendix_c.json` |
| Fresh execution of Section 4 notebooks | `sage -python certificates/notebooks.py` | executed notebooks and `results/recovered_parameters.json` |
| Recurrences and first 20 coefficients | `sage -python certificates/coefficients.py --count 20 --verify-count 32` | `results/coefficients.json` |
| Exact quaternion and CM input validation | `sage -python certificates/nelson_inputs.py` | `results/nelson_inputs.json` |
| Nelson values and cutoff comparison | `sage -python certificates/numerical.py` | `results/numerical.json` |
| Cross-certificate consistency, after notebooks and Appendix C | `sage -python certificates/integration.py` | `results/integration.json` |
| Four displayed CM-image polynomial sanity checks | `sage -python certificates/cm_claims.py` | `results/cm_claims.json` |
| Twelve-under-200 setup sanity only | `sage certificates/cm_sanity.sage` | `results/cm_sanity.json` |

## Scope and provenance

`inputs/` contains snapshots of the supplied draft and notebooks. The
executed notebooks are separate files in `results/`. The Nelson evaluator
is accompanied by exact-input and numerical-comparison wrappers. Published
coordinate values and theorem inputs are recorded in
[data/appendix_c_sources.json](data/appendix_c_sources.json).

The under-200 directory retains the supplied reproduction code and exact
CM input data. Two API compatibility adaptations make its sanity check run
under Sage 9.5: the fundamental-discriminant predicate is expressed as
`fundamental_discriminant(D) == D`, and quaternion inversion uses `~q`.
The full batch, analytic CM-image integration, and image recognition were
not rerun. The four displayed Table 4 polynomials have independent exact
sanity checks; see [their scope](docs/cm_claims.md).

The original external files and the draft text were not edited. This package
contains the directly relevant source, data, documentation, and execution
records; previous exploratory work is not required.
