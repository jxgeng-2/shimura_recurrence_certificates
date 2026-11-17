# Table 3 certificate

Run from the package directory:

```console
python3 certificates/families.py
```

The script uses exact SymPy rational arithmetic. Its report is
`data/families_certificate.json`. All twelve rows pass.

## Inputs and logical scope

The inputs are the nine algebraic models and eigendifferentials in Table 1,
the quotient coordinates, and the elliptic fibers identified in Appendix C.
The script does not infer elliptic orders from a plane curve: those are
uniformization data and remain explicit geometric inputs. Their coordinate
identifications belong to the Appendix C certificate. No recognized numerical
Taylor coefficient is an input to this computation.

The family and H columns of the draft are separately transcribed as expected
answers. They are compared with the derived answers only after the derivation.

## What is derived and checked

For each form with model `y^2 = F(x)` and differential `omega = h(x) dx/y`,
the script computes

\[
H(u(x))=\frac{u_x(x)^2 F(x)}{h(x)^2}.
\]

It verifies this identity as an exact rational function. For level 34 it
descends the computed rational function using `x^2 - u*x - 1 = 0`;
it does not assume the displayed formula for H. The model is checked to be
squarefree, the involutions are checked to preserve it and the quotient
coordinate, and the differential is checked to be an eigenvector. The
degree-two map from the x-line gives the degree-four quotient map from the
hyperelliptic curve.

The prescribed polar coefficient is computed by differentiating the local
models `z=v^(1/e)` and `t=alpha*v^(m+1)/(m+1)`. The Schwarzian chain rule gives

\[
\{z,t\}=\frac{\{z,v\}-\{t,v\}}{(dt/dv)^2}
=\frac{(m+1)^2-e^{-2}}{2\alpha^2}v^{-2m-2}.
\]

This highest polar term is unchanged by the higher analytic local terms.
All four combinations `(m,e)` required by the table are checked directly.
The script then applies the local coordinate change to obtain:

| Fiber type | Prescribed coefficient |
| --- | --- |
| Unramified elliptic, `(m,r)=(0,1)` | `(u-a)^(-2)`: `kappa_e H(a)` |
| Ramified elliptic, `(m,r)=(0,2)` | `(u-a)^(-1)`: `kappa_e H'(a)/4` |
| Differential zero, `(m,r)=(1,2)` | `(u-a)^(-2)`: `kappa_(2e) H(a)` |
| Infinity, `(m,r)=(0,2)` | `u`: `kappa_e a_N` |

Here `kappa_e=(1-e^(-2))/2`. At every zero fiber the script checks
`u=v^2`, the simple zero of the differential, the absence of its next odd
correction, and `H(0)=4/alpha^2`. For infinity it computes `(dv/dt)^2`
from the actual model and verifies its nonzero finite limit. These local
checks also verify the asserted ramification types.

The free coefficients are obtained from

\[
\Delta_\omega=\sum_a
\left\lfloor\frac{2m_a+\mathbf 1_{e_a>1}}{r_a}\right\rfloor[a].
\]

The basis is `1` together with `1/(u-a)` for each point of this divisor.
The script verifies its linear independence exactly and checks its dimension
against the last column of Table 3. Adding its arbitrary linear combination
to the prescribed principal parts produces the full family; the resulting
identity with the draft is checked for each of the twelve forms. The fixed
principal parts are also rechecked by exact limits.

The dimensions in table order are

```text
14a 15a 21a 33a 34a 46a 26a 26b 38a 38b 58a 58b
 1   1   2   2   2   2   2   2   2   2   3   3
```

## Reuse

`derive_H(form)` returns the derived H and its model metadata.
`derive_family(form)` returns one complete certificate row.
`run_certificate(output_dir=None)` returns the complete JSON-serializable
report; if an output directory is supplied, it writes
`families_certificate.json` there.

Every mathematical check uses an explicit exception on failure and remains
active under `python -O`.
