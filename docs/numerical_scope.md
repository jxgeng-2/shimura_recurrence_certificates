# What the numerical certificate establishes

`certificates/numerical.py` evaluates every assumed initial coefficient and
every CM expansion scale for the twelve forms. It records the largest observed
differences in `results/numerical.json`, separately for:

| Quantity | Report field |
| --- | --- |
| Normalized coefficient versus the stated rational value | `max_coefficient_discrepancy` |
| Theta versus its stated period times algebraic value | `max_theta_discrepancy` |
| Theta/(pi Omega_D^2) versus its stated algebraic value | `max_normalized_alpha_discrepancy` |
| Change when the lattice cutoff is increased | `max_cutoff_difference` |

The numerical checks require differences below `10^-40`, as requested. The
recorded values give the actual achieved agreement, rather than just that
pass/fail threshold. The draft uses `10^-45` for initial coefficients and
`10^-50` for the normalized scale; these are different comparisons and should
be read from their respective report fields. The completed report contains
the individual values, imaginary parts, cutoffs, and comparisons for each form.

The standard run uses 90 decimal digits, an 80-digit target for the weight
function quadrature, and squared lattice cutoffs 144 and 196. It also compares
the eta-product period with the Gamma-product period at high precision.

## Exact statements and numerical statements

The quaternion multiplication, lattice bases, trace pairings, reduced order
discriminants, CM discriminants, and grouped lattice weights use exact
integer or rational arithmetic. The additional checks verify that the orders
contain 1 and are closed under multiplication. Each stated tau is checked
against the quadratic fixed-point equation in an exact radical algebra and
against the displayed expression in the appendix.

`certificates/nelson_inputs.py`, run with `sage -python`, independently checks
the algebra discriminant, maximality of each order, the trace-dual bases,
matrix multiplication in the real embedding, the maximal CM order and its
class number one, and the exact positivity of the imaginary part of tau.
It also checks all twelve elliptic conductors and their first 100 Hecke
coefficients against Sage. All these checks pass; their complete record is
`results/nelson_inputs.json`. This script evaluates no numerical lattice sums.

The evaluations of the infinite Nelson sum are numerical. The program sums
a finite set of lattice points, approximates the weight-function integrals,
and performs rounded decimal arithmetic. It reports observed agreement and
stability under an increased cutoff. It does not attach an interval enclosure
to the omitted infinite tail, quadrature remainder, or propagated rounding
through the ratios involving the zeroth and first sums. Consequently the
reported numerical differences do not supply a rigorous upper bound for the
distance between the analytic coefficient and the conjectured exact value.
The recurrence calculations retain the draft's conditional interpretation.

The implemented specialization is consistent with Nelson's coefficient
factor, locally dual lattice, matrix coordinates, and derivative weight for
weight two and derivative indices `(0,n)`. In the present trivial-level
setting the coefficient factor is `1/(d*a_f(d))`. Nelson's convergence
argument supports truncation of the sum; an explicit numerical enclosure
requires additional quantitative estimates. See
[Nelson, Theorems 3.1 and 3.9 and Remark 3.3](https://arxiv.org/pdf/1210.1243).

## Interpreting tau and the scale

The exact tau check identifies the chosen fixed point of the specified CM
embedding. The quaternion presentations and maximal orders are explicit
normalization data. At level 34, the supplied order was chosen to match the
draft's scale convention; it is recorded as part of that convention. Numerical
agreement with the scale is evidence for the stated normalization, not an
independent proof of the algebraic scale formula.

For a future rigorous analytic error bound, the additional ingredients are
an explicit lattice-tail estimate, enclosed evaluations of the weight
integrals and periods, and nonzero enclosures for the denominator sums before
forming coefficient ratios. None of these is needed for the exact algebraic
family and recurrence certificates.
