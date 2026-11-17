# Exact CM expansion coefficients

`certificates/coefficients.py` computes the first N coefficients for each of
14a, 15a, 21a, 33a, 34a, 46a, 26a, 26b, 38a, 38b, 58a, and 58b using only
Python's standard library. All arithmetic is exact rational arithmetic.

From the package directory:

```python
from certificates.coefficients import first_coefficients, appendix_expansion, w_coefficients

c = first_coefficients("14a", 20)  # exactly c_0,...,c_19, as Fraction objects
expansion = appendix_expansion("14a", 20)
full_w = w_coefficients("14a", 20) # coefficients of f, including (1-w)^2
```

The appendix convention is

\[
w=\frac{z-\tau}{z-\bar\tau},\quad q=\Theta w,\qquad
f(z)=(1-w)^2 A(\Theta w),\qquad A(q)=\sum_{n\ge0}c_n\frac{q^n}{n!}.
\]

Thus the returned `c` list is the list printed before the factorials in the
appendix; the ordinary coefficients of A are `c[n]/factorial(n)`.
`appendix_expansion` also supplies the draft's expressions for tau and
Theta/(pi Omega_D^2). These expressions are metadata, not a new proof of their
arithmetic identification. Principal square roots are intended, and the
displayed tau is the upper-half-plane representative.

`w_coefficients(form, N)` supplies polynomials in the formal variable Theta,
represented as dictionaries from degree to exact rational coefficient.
Alternatively, supply a numerical or exact `theta` argument to evaluate them.
If `a_n=c_n/n!`, the coefficient of w^n in f is

\[
a_n\Theta^n-2a_{n-1}\Theta^{n-1}+a_{n-2}\Theta^{n-2},
\]

with negative subscripts interpreted as zero. This explicitly includes the
prefactor `(1-w)^2`; it does not confuse the coefficients of f with those of A.

The recurrence uses the two right-hand sides

\[
A''=\frac{3(A')^2}{2A}-L R(u)A^3,\qquad
u''=\frac{A'}A u'+\frac L2 A^2H'(u).
\]

At step k, their q^k coefficients depend only on known coefficients through
k+1. Division by `(k+1)(k+2)` gives the next ordinary coefficients. All
denominators evaluated at the twelve centers are nonzero, and A(0)=1, so
the process determines a unique formal solution over Q for every requested
finite N. Arbitrarily large N is permitted; no asymptotic performance claim
is made for this straightforward exact series implementation.

Run the certificate with:

```text
python certificates/coefficients.py --count 20 --verify-count 32
```

It prints all first 20 coefficients and writes `results/coefficients.json`
and `results/coefficients_first20.md`. It computes 32 coefficients for each
form and verifies the following exact identities through q^29:

\[
AA''-\tfrac32(A')^2+LR(u)A^4=0,\qquad
Au''-A'u'-\tfrac L2 A^3H'(u)=0,\qquad
(u')^2-LA^2H(u)=0.
\]

The last identity checks the first integral, including the initial branch
condition u_1^2=L H(u_0), independently of the divided update rule. The
certificate also compares every coefficient stated in Appendices A and B
with the generated value and verifies multiplication by `(1-w)^2` at an
exact test value of Theta. The 32-coefficient verification is an exact
finite calculation; the recursive construction explains extension to any N.

These are certificates for the conditional theorem. They do not turn the
numerically recognized initial coefficients or periods into unconditional
algebraic results. A separate numerical computation is needed to support
those input recognitions.
