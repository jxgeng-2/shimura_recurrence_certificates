# Appendix C: exact CM-coordinate certificate

This certificate covers **all 24 rows of Table 7** in the supplied draft, the three marked reciprocal conventions, and the three-center argument in (C.1). Run `python certificates/appendix_c.py`. The output is `results/appendix_c.json`; no floating-point arithmetic or numerical recognition is used.

## Published coordinates and explicit changes

The source locations are recorded in `data/appendix_c_sources.json`. Elkies uses coordinates `(s,v)` with `T=-s^2` at level 14 and `T=-3s^2` at level 15. The full model changes are `(x,y)=(4s,4v)` at level 14 and `(x,y)=(s,v)` at level 15. Substitution gives the draft's quartics identically, hence `u=-16T` and `u=-T/3`.

For Lin--Yang's marked coordinate `j`, let `F_N(j)` denote the branch polynomial read from their Table 7. The program verifies each of the following identities with the draft's quartic `P_N`:

| N | x(j) | F_N(j) | P_N(x(j))/F_N(j) |
|---|---|---|---|
| 21 | 7/j | 7j^4-94j^2+343 | -49/j^4 |
| 33 | 1/j | 243j^4+10j^2+3 | -1/j^4 |
| 34 | (5+3j)/(3-5j) | 27j^4+136j^3+122j^2-136j+27 | -289/(3-5j)^4 |
| 46 | 8j | 64j^4-45j^2+8 | -64 |

These changes lift explicitly from the auxiliary equation `Y^2=-F_N(j)` by `y=7Y/j^2`, `Y/j^2`, `17Y/(3-5j)^2`, and `8Y`, respectively. The branch-polynomial computation verifies the coordinate algebra; identifying the marked quotient with the Shimura curve uses the cited model theorem and normalizations. A branch divisor alone does not determine the arithmetic twist.

The quotient coordinates are `49/j^2`, `1/j^2`, `x(j)-1/x(j)`, and `64j^2`. The program verifies invariance under the published Atkin--Lehner actions `j -> -j`, except at level 34 where the action is `j -> -1/j`. It evaluates the published marked values at finite points and at infinity. Thus reciprocal or sign choices cannot be silently changed.

For levels 26 and 38, Guo--Yang use `s=x^2`; for level 58 their Table 16 uses **`s=-x^2/19`**. The program applies these changes before taking `(u_a,u_b)=(x^2,x^-2)`. In particular `s=11/19` and `27/19` at level 58 become `(-11,-1/11)` and `(-27,-1/27)`. It includes all zero and infinite fibers. The published CM identifications are explicit inputs, not asserted as a new implementation of Borcherds products.

## The Lin--Yang intersection calculation

The script includes all four complete polynomials from Lin--Yang Table 4 for each of levels 21, 33 and 46. It substitutes the draft's rational `j` values and checks weighted projective equality using weights `(2,3,5,6)`. In ordinary coordinates the exact calculation is

`S_N(j)_k = lambda^k T_k`,

where `T=(a^2,b^2,a^2 b^2,a^6+b^4)` and:

| N | D | j | a | b | lambda |
|---|---|---|---|---|---|
| 21 | -43 | 5/3 | 10 | 21 | 16/9 |
| 33 | -67 | 5/3 | 20 | 231 | 16/9 |
| 46 | -163 | 7/9 | 1870 | 27531 | 16/19683 |

Writing `q=b^4/a^6` and `mu=a^4/b^2`, the separate identity `T_k=mu^k S_6(q)_k` is checked in each of the four coordinates. All `s5` values are nonzero, so these points lie off the decomposable locus. The complete rational coordinates and scales are included in the machine-readable results.

Errthum's normalization and Lin--Yang's three marked values give `J=-16t6/27`. The certificate independently transcribes both norm factorizations in Errthum Table 3 and checks `|t6|=27q/16` and `|1-t6|=1+27q/16`. If `A=27q/16`, the two equations give

`Re(t6)=(1+A^2-(1+A)^2)/2=-A`,

`Im(t6)^2=A^2-Re(t6)^2=0`.

Thus `t6=-27q/16` and `J=q`, even without a preliminary assumption that `t6` is real. The positive signs shown for these entries in Table 2 of the linked Errthum version are not used; Table 3's pair of norms resolves the sign.

## Why the common surface has the asserted CM discriminant

This step is a mathematical argument using the moduli interpretation, rather than a numerical test. The weighted invariant identity identifies the same principally polarized complex abelian surface in the discriminant-6 and discriminant-N moduli images. On the discriminant-6 side, the published norm theorem identifies its CM order as the maximal order `O_D`.

At a QM CM point the commuting `O_D` action on the tangent space is scalar. Its period lattice is therefore a torsion-free rank-two `O_D` module. For `D=-43,-67,-163` the script enumerates all reduced primitive positive definite forms, obtaining respectively `(1,1,11)`, `(1,1,17)`, `(1,1,41)`. Thus the class number is one. The module is free, so, after a complex linear change of basis, the unpolarized surface is `E_D^2` and its endomorphism ring is `M_2(O_D)`.

The other quaternion algebra `B_N` acts inside `M_2(K)`, where `K=Q(sqrt(D))`. Its centralizer is the scalar `K`: extending scalars to `K` identifies the four-dimensional algebra `B_N tensor K` with `M_2(K)`. Intersecting the centralizer with `M_2(O_D)` gives exactly `O_D`, not an unspecified suborder. Hence the point on the discriminant-N curve has CM by the same maximal order and therefore discriminant `D`. The script also checks the nonsplitting conditions at every prime dividing N.

Finally Gonzalez--Rotger Proposition 5.6 makes each of these class-number-one CM sets a single Atkin--Lehner orbit. The quotient coordinates are invariant on that orbit, and the verified substitutions give `441/25`, `9/25`, and `3136/81`.

This argument is confined to the three selected centers. The extra table discriminants `-52` and `-232` have class number **two**, as the certificate also checks. Their zero-fiber identifications come from Guo--Yang, not from a false class-number-one assertion.

## Source boundary

The executable part proves the rational algebra, weighted invariant equalities, norm-sign recovery, and reduced-form counts. The moduli interpretation, published CM values, and transitivity theorem remain cited mathematical inputs. This is an exact certificate relative to those theorems, not a formal proof assistant reconstruction of the literature. It does not claim that substituting into a curve equation alone identifies a point as CM.

Sources: [Elkies](https://arxiv.org/pdf/math/0005160), [Lin--Yang](https://arxiv.org/pdf/1807.00466), [Guo--Yang](https://arxiv.org/pdf/1510.06193), [Errthum](https://arxiv.org/pdf/0711.4316), [Gonzalez--Rotger](https://arxiv.org/pdf/math/0612732). The table/example numbers refer to the expanded versions linked here.
