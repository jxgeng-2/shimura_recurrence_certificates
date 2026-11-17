#!/usr/bin/env python3
"""
nelson_cm_values.py -- plain numerical evaluation (no certified error bounds) of the CM data of
"Recurrence formulas for CM expansions on genus-one Shimura curves" by Nelson's formula:

  * the normalized coefficients c_n of Table 4, aiming at |c_num - c_table| < 1e-40;
  * the exact CM point tau_N of Table 5 (closed form, computed from the quaternion data);
  * Theta_N = b_1/b_0, compared with the algebraic value pi Omega_D^2 * alpha_N of Table 5 to < 1e-40.

    python3 nelson_cm_values.py            # all twelve forms   (about 3 minutes)
    python3 nelson_cm_values.py 14 46      # only these levels
    python3 nelson_cm_values.py --orders 34   # list the maximal orders containing Z<i, j, mu> and their Theta

Needs Python 3 and numpy (point counts).  Uses exact integers/fractions for every lattice sum and the
decimal module (70 digits) for the transcendental part.

Conventions of Table 5.  B = (a, b)_Q with i^2 = a, j^2 = b, k = ij.  For a > 0 > b,
    x + y i + z j + w k  ->  [[x + y sqrt a, (z + w sqrt a) sqrt|b|], [-(z - w sqrt a) sqrt|b|, x - y sqrt a]],
and for a < 0 < b the same with i and j interchanged (k -> -k):
    x + y i + z j + w k  ->  [[x + z sqrt b, (y - w sqrt b) sqrt|a|], [-(y + w sqrt b) sqrt|a|, x - z sqrt b]].
mu in R has trace 1 and discriminant D; tau_N is its fixed point in H.  The data (a, b, R, mu) below reproduce every
tau_N of Table 5 exactly.  N = 26, 38: Magma's maximal order (generators of Gamma) and mu from the user's notebook;
N = 58: the order constructed in the earlier genus-two session.  N = 14, 15, 21, 33, 34, 46: (a, b, mu) are read off
from the printed tau_N, and R is a maximal order containing Z<i, j, mu> (Z-bases below) -- the only one, except for
N = 34 (two such orders), where the one that reproduces Theta_34 of Table 5 is taken (see --orders).  The c_n do not depend on any of these choices;
Theta_N does (through a factor of absolute value 1), which is why the data must match the table's conventions.

Formula.  P. D. Nelson, Math. Comp. 84 (2015), Thms 3.1 and 3.9 with z1 = z2 = tau, (j1, j2) = (0, n), k = 2:
    y^{2+n} conj(f(tau)) delta^n f(tau) = sum_{d|N} c_f(d) sum_{beta in R^(d), nrd beta > 0} a_f(d nrd beta)
                                        (-4 pi)^n (X/|X|)^{2+n} Y^n W_{2+n}(4 pi |X|),
    c_f(d) = 1/(d a_f(d)),   R^(d) = R + (N/d) R^*,   W_k(x) = 2^{1-k} sum_m m^k (m x K_{k-1}(m x) - K_k(m x)),
with X, Y Nelson's coordinates of sigma_tau^{-1} beta sigma_tau.  Writing beta = kappa1 + kappa2 j0 with
kappa_i in K = Q(mu) and j0 in K^perp:  X = iota(kappa1),  Y = conj(iota(kappa2)) Y0,  Y0 = Y(j0),  where
iota: K -> C is the embedding fixing tau.  With S'_n the sum without the factor (-4 pi Y0)^n:
    b_n/b_0 = (-4 pi Y0)^n S'_n / (n! S'_0),   Theta = b_1/b_0 = -4 pi Y0 S'_1/S'_0,   c_n = S'_n S'_0^{n-1} / S'_1^n.
Periods: pi Omega_D^2 with Omega_D^2 = |eta((1+sqrt D)/2)|^4 (Chowla-Selberg, h(D) = 1), checked against the
Gamma-product definition of the paper in double precision.
"""
import sys
import math
import time
import itertools
from fractions import Fraction as Fr
from functools import reduce
from collections import defaultdict
from decimal import Decimal, getcontext
import numpy as np

DIGITS = 70                 # working precision (decimal digits)
MCUT = 144                  # lattice cutoff |X|^2 <= MCUT
WDIG = 60                   # target relative accuracy of W_k
getcontext().prec = DIGITS


# =========================================================================================== quaternion algebra
class QuatAlg:
    def __init__(self, a, b):
        self.a, self.b = Fr(a), Fr(b)

    def mul(self, x, y):
        a, b = self.a, self.b
        x1, y1, z1, w1 = x
        x2, y2, z2, w2 = y
        return (x1 * x2 + a * y1 * y2 + b * z1 * z2 - a * b * w1 * w2,
                x1 * y2 + y1 * x2 - b * z1 * w2 + b * w1 * z2,
                x1 * z2 + z1 * x2 + a * y1 * w2 - a * w1 * y2,
                x1 * w2 + w1 * x2 + y1 * z2 - z1 * y2)

    def nrd(self, x):
        a, b = self.a, self.b
        return x[0] ** 2 - a * x[1] ** 2 - b * x[2] ** 2 + a * b * x[3] ** 2

    def pair(self, x, y):                      # trd(x conj(y))
        a, b = self.a, self.b
        return 2 * (x[0] * y[0] - a * x[1] * y[1] - b * x[2] * y[2] + a * b * x[3] * y[3])


def lcm(a, b):
    return a * b // math.gcd(a, b)


def hnf_rows(M):
    """Row Hermite normal form (upper triangular, positive pivots); nonzero rows."""
    A = [list(r) for r in M]
    m, n = len(A), len(A[0])
    r = 0
    for c in range(n):
        while True:
            nz = [i for i in range(r, m) if A[i][c] != 0]
            if not nz:
                break
            piv = min(nz, key=lambda i: abs(A[i][c]))
            A[r], A[piv] = A[piv], A[r]
            done = True
            for i in range(r + 1, m):
                if A[i][c]:
                    q = A[i][c] // A[r][c]
                    A[i] = [x - q * y for x, y in zip(A[i], A[r])]
                    done = done and A[i][c] == 0
            if done:
                break
        if r < m and A[r][c] != 0:
            if A[r][c] < 0:
                A[r] = [-x for x in A[r]]
            for i in range(r):
                q = A[i][c] // A[r][c]
                A[i] = [x - q * y for x, y in zip(A[i], A[r])]
            r += 1
            if r == m:
                break
    return [row for row in A if any(row)]


def lattice_basis(vecs):
    vecs = [tuple(Fr(x) for x in v) for v in vecs]
    den = reduce(lcm, [x.denominator for v in vecs for x in v], 1)
    return [tuple(Fr(x, den) for x in row) for row in hnf_rows([[int(x * den) for x in v] for v in vecs])]


def det(M):
    M = [[Fr(x) for x in r] for r in M]
    n, d = len(M), Fr(1)
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d


def inv(M):
    n = len(M)
    A = [[Fr(x) for x in r] + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        p = A[c][c]
        A[c] = [x / p for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [r[n:] for r in A]


def contains(basis, v):
    Bi = inv([list(b) for b in basis])
    return all(sum(v[k] * Bi[k][j] for k in range(4)).denominator == 1 for j in range(4))


def ring_closure(A, gens, maxit=20):
    basis = lattice_basis([(1, 0, 0, 0)] + [tuple(Fr(x) for x in g) for g in gens])
    for _ in range(maxit):
        nb = lattice_basis(basis + [A.mul(x, y) for x in basis for y in basis])
        if nb == basis:
            return basis
        basis = nb
    raise RuntimeError("closure did not stabilise")


def reduced_discriminant(A, basis):
    d = abs(det([[A.pair(x, y) for y in basis] for x in basis]))
    rn, rd = math.isqrt(d.numerator), math.isqrt(d.denominator)
    assert rn * rn == d.numerator and rd * rd == d.denominator
    return Fr(rn, rd)


def dual(A, basis):
    Gi = inv([[A.pair(x, y) for y in basis] for x in basis])
    return [tuple(sum(Gi[i][l] * basis[l][c] for l in range(4)) for c in range(4)) for i in range(4)]


def maximal_orders_containing(A, gens, N):
    """All maximal orders (reduced discriminant N) containing the ring generated by gens (p-enlargements)."""
    def integral(x):
        return (2 * x[0]).denominator == 1 and A.nrd(x).denominator == 1
    todo, found, seen = [ring_closure(A, gens)], [], set()
    while todo:
        O = todo.pop()
        if tuple(O) in seen:
            continue
        seen.add(tuple(O))
        disc = reduced_discriminant(A, O)
        if disc == N:
            if not any(all(contains(U, v) for v in O) for U in found):
                found.append(O)
            continue
        m = int(disc / N)
        p = next(q for q in range(2, m + 1) if m % q == 0)
        for c in itertools.product(range(p), repeat=4):
            if not any(c):
                continue
            x = tuple(sum(Fr(ci, p) * O[k][t] for k, ci in enumerate(c)) for t in range(4))
            if not integral(x):
                continue
            try:
                O2 = ring_closure(A, list(O) + [x], maxit=6)
            except RuntimeError:
                continue
            if all(integral(v) for v in O2) and all(integral(A.mul(u, v)) for u in O2 for v in O2):
                d2 = reduced_discriminant(A, O2)
                if d2 < disc and (d2 / N).denominator == 1:
                    todo.append(O2)
    return found


# =========================================================================================== a_n of elliptic curves
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0].tolist()


def ap_curve(ainv, p):
    """a_p = p + 1 - #E(F_p) counting all projective points (also correct, = 1/-1/0, at bad p)."""
    a1, a2, a3, a4, a6 = ainv
    if p == 2:
        cnt = 1 + sum(1 for x in range(2) for y in range(2)
                      if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % 2 == 0)
        return p + 1 - cnt
    b2, b4, b6 = a1 * a1 + 4 * a2, a1 * a3 + 2 * a4, a3 * a3 + 4 * a6
    x = np.arange(p, dtype=np.int64)
    rhs = (((4 * x % p) * x % p) * x + (b2 % p) * (x * x % p) + (2 * b4 % p) * x + b6) % p
    sq = np.zeros(p, dtype=np.int64)
    sq[(x * x) % p] = 1
    return -int(np.where(rhs == 0, 0, np.where(sq[rhs] == 1, 1, -1)).sum())


def an_table(ainv, nmax, bad):
    a = [0] * (nmax + 1)
    a[1] = 1
    ps = primes_upto(nmax)
    pp = {}
    for p in ps:
        ap = ap_curve(ainv, p)
        vals, q, k = {0: 1, 1: ap}, p, 1
        while q * p <= nmax:
            q, k = q * p, k + 1
            vals[k] = vals[k - 1] * ap if p in bad else ap * vals[k - 1] - p * vals[k - 2]
        pp[p] = vals
    spf = list(range(nmax + 1))
    for p in ps:
        if p * p > nmax:
            break
        for m in range(p * p, nmax + 1, p):
            if spf[m] == m:
                spf[m] = p
    for n in range(2, nmax + 1):
        p, m, k = spf[n], n, 0
        while m % p == 0:
            m, k = m // p, k + 1
        a[n] = pp[p][k] * a[m]
    return a


# =========================================================================================== exact radical expressions
def squarefree(m):
    """m = s^2 * t with t squarefree (sign kept in t)."""
    s, t, f = 1, abs(m), 2
    while f * f <= t:
        while t % (f * f) == 0:
            t //= f * f
            s *= f
        f += 1
    return s, (t if m > 0 else -t)


class Rad:
    """Finite sums  sum c_m sqrt(m)  (c_m rational, m squarefree, m < 0 meaning i sqrt|m|)."""

    def __init__(self, d=None):
        self.d = {}
        for m, c in (d or {}).items():
            s, t = squarefree(m)
            if c:
                self.d[t] = self.d.get(t, Fr(0)) + Fr(c) * s
        self.d = {m: c for m, c in self.d.items() if c}

    def __add__(self, o):
        r = dict(self.d)
        for m, c in o.d.items():
            r[m] = r.get(m, Fr(0)) + c
        return Rad(r)

    def __mul__(self, o):
        if not isinstance(o, Rad):
            return Rad({m: c * Fr(o) for m, c in self.d.items()})
        r = Rad()
        for m1, c1 in self.d.items():
            for m2, c2 in o.d.items():
                sign = -1 if (m1 < 0 and m2 < 0) else 1
                r = r + Rad({m1 * m2 if sign == 1 else m1 * m2: sign * c1 * c2})
        return r

    def __eq__(self, o):
        return self.d == o.d

    def value(self):
        """(re, im) as Decimals."""
        re, im = Decimal(0), Decimal(0)
        for m, c in self.d.items():
            v = Decimal(c.numerator) / Decimal(c.denominator) * Decimal(abs(m)).sqrt()
            if m > 0:
                re += v
            else:
                im += v
        return re, im

    def __str__(self):
        if not self.d:
            return '0'
        den = reduce(lcm, [c.denominator for c in self.d.values()], 1)
        out = ''
        for m in sorted(self.d, key=lambda t: (t < 0, abs(t))):
            num = int(self.d[m] * den)
            mag = '' if (abs(num) == 1 and m != 1) else str(abs(num))
            rad = '' if m == 1 else ('sqrt(%d)' % m)
            body = mag + ('*' if mag and rad else '') + rad
            out += ('-' if num < 0 else '') + body if not out else (' - ' if num < 0 else ' + ') + body
        return '(%s)/%d' % (out, den) if den != 1 else out


# =========================================================================================== data (Tables 1, 4, 5)
CURVES = {'14a': [1, 0, 1, 4, -6], '15a': [1, 1, 1, -10, -10], '21a': [1, 0, 0, -4, -1], '33a': [1, 1, 0, -11, 0],
          '34a': [1, 0, 0, -3, 1], '46a': [1, -1, 0, -10, -12], '26a': [1, 0, 1, -5, -8], '26b': [1, -1, 1, -3, 3],
          '38a': [1, 0, 1, 9, 90], '38b': [1, 1, 1, 0, 1], '58a': [1, -1, 0, -1, 1], '58b': [1, 1, 1, 5, 9]}

h = Fr(1, 2)
DATA = {   # a, b, generators of R (as a ring; a Z-basis is fine too), mu, D, forms
    14: dict(a=7, b=-2, gens=[(h, h, 0, h), (0, 1, 0, 0), (0, 0, h, h), (0, 0, 0, 1)],
             mu=(h, -h, Fr(-3, 2), 0), D=-11, forms=['14a']),
    15: dict(a=15, b=-30, gens=[(h, 0, Fr(1, 4), Fr(1, 20)), (0, h, Fr(1, 4), Fr(1, 60)), (0, 0, h, Fr(1, 30)),
                                (0, 0, 0, Fr(1, 15))],
             mu=(h, 0, Fr(-1, 4), Fr(1, 60)), D=-7, forms=['15a']),
    21: dict(a=21, b=-1, gens=[(h, h, 0, 0), (0, 1, 0, 0), (0, 0, h, h), (0, 0, 0, 1)],
             mu=(h, -h, -10, -2), D=-43, forms=['21a']),
    33: dict(a=33, b=-1, gens=[(h, h, 0, 0), (0, 1, 0, 0), (0, 0, h, h), (0, 0, 0, 1)],
             mu=(h, -h, -5, 0), D=-67, forms=['33a']),
    34: dict(a=34, b=-3, gens=[(h, 0, Fr(1, 6), Fr(1, 3)), (0, h, 0, h), (0, 0, Fr(1, 3), Fr(2, 3)), (0, 0, 0, 1)],
             mu=(h, h, Fr(-7, 2), h), D=-11, forms=['34a']),
    46: dict(a=46, b=-3, gens=[(h, 0, Fr(1, 6), Fr(1, 3)), (0, h, 0, h), (0, 0, Fr(1, 3), Fr(2, 3)), (0, 0, 0, 1)],
             mu=(h, 0, Fr(-35, 6), Fr(-2, 3)), D=-163, forms=['46a']),
    26: dict(a=26, b=-5, gens=[(5, -1, Fr(6, 5), Fr(1, 5)), (5, -1, Fr(-6, 5), Fr(-1, 5)), (Fr(5, 2), -h, -h, 0),
                               (Fr(3, 2), 0, Fr(1, 10), Fr(1, 10)), (5, 1, Fr(-6, 5), Fr(-1, 5)), (5, 1, Fr(6, 5), Fr(1, 5)),
                               (Fr(5, 2), h, h, 0), (Fr(3, 2), 0, Fr(-1, 10), Fr(-1, 10))],
             mu=(h, 1, Fr(-7, 2), -h), D=-11, forms=['26a', '26b']),
    38: dict(a=-19, b=2, gens=[(3, 0, -2, 0), (2, -1, Fr(-5, 2), h), (Fr(5, 2), h, -h, -h), (0, -1, -3, 0),
                               (Fr(5, 2), -h, -h, -h), (3, 0, 2, 0), (2, 1, Fr(5, 2), -h), (Fr(5, 2), -h, h, h),
                               (0, 1, 3, 0), (Fr(5, 2), h, h, h)],
             mu=(h, -h, -1, 0), D=-11, forms=['38a', '38b']),
    58: dict(a=-29, b=2, gens=[(h, h, 0, h), (0, 1, 0, 0), (0, 0, h, h), (0, 0, 0, 1)],
             mu=(h, h, Fr(3, 2), 0), D=-11, forms=['58a', '58b']),
}

TABLE4 = {
    '14a': {2: Fr(-4), 3: Fr(214, 9)},
    '15a': {2: Fr(-1, 4), 3: Fr(27, 16)},
    '21a': {2: Fr(-905, 3969), 3: Fr(-1318, 27783)},
    '33a': {2: Fr(-72337, 529), 3: Fr(9202302, 12167), 4: Fr(3009372468, 279841), 5: Fr(-1149714703920, 6436343)},
    '34a': {2: Fr(-12, 25), 3: Fr(2762, 125), 4: Fr(-17736, 625), 5: Fr(-67648, 125)},
    '46a': {2: Fr(27824049, 123402050), 3: Fr(-114816966372, 484661551375),
            4: Fr(-784666539922896, 3807016486050625), 5: Fr(-3694325348311426944, 5980822899585531875)},
    '26a': {2: Fr(-16, 147), 3: Fr(274, 3087), 4: Fr(-2488, 9261), 5: Fr(352768, 1361367)},
    '26b': {2: Fr(40), 3: Fr(-734), 4: Fr(-11240), 5: Fr(1308544)},
    '38a': {2: Fr(-64, 27), 3: Fr(-62, 27), 4: Fr(-19672, 2187), 5: Fr(41984, 6561)},
    '38b': {2: Fr(-16, 169), 3: Fr(-878, 2197), 4: Fr(7688, 28561), 5: Fr(-2699648, 371293)},
    '58a': {2: Fr(-4, 5), 3: Fr(314, 125), 4: Fr(31192, 625), 5: Fr(-1977728, 625), 6: Fr(4428752, 3125)},
    '58b': {2: Fr(140, 9), 3: Fr(262, 3), 4: Fr(-39592, 81), 5: Fr(420224, 81)},
}

TABLE5_TAU = {
    14: Rad({14: Fr(-1, 6), -22: Fr(1, 6)}),
    15: Rad({-210: Fr(1, 14), -14: Fr(-1, 14)}),
    21: Rad({1: Fr(-21, 16), 21: Fr(-5, 16), -43: Fr(5, 16), -903: Fr(1, 16)}),
    33: Rad({33: Fr(-1, 10), -67: Fr(1, 10)}),
    34: Rad({3: Fr(-34, 45), 102: Fr(7, 45), -33: Fr(7, 45), -1122: Fr(-1, 45)}),
    46: Rad({-489: Fr(35, 489), -22494: Fr(4, 489)}),
    26: Rad({5: Fr(52, 115), 130: Fr(14, 115), -55: Fr(7, 115), -1430: Fr(1, 115)}),
    38: Rad({38: Fr(-2, 19), -209: Fr(1, 19)}),
    58: Rad({58: Fr(-3, 29), -319: Fr(1, 29)}),
}

_X26 = Rad({1: 988, 26: -33, -11: -156, -286: -19})
TABLE5_THETA = {   # Theta_N / (pi Omega_D^2)
    '14a': Rad({1: 1, -77: Fr(1, 7)}),
    '15a': Rad({15: Fr(-8, 15)}),
    '21a': Rad({1: Fr(-63, 4), -903: Fr(-3, 4)}),
    '33a': Rad({1: Fr(-23, 5), -2211: Fr(-23, 165)}),
    '34a': Rad({1: 1, 34: Fr(-11, 102), -11: Fr(-1, 3), -374: Fr(-1, 34)}),
    '46a': Rad({46: Fr(7855, 69)}),
    '26a': _X26 * Fr(-21, 2990),
    '26b': _X26 * Fr(1, 2990),
    '38a': Rad({1: 4, -22: 1}) * Fr(-9, 19),
    '38b': Rad({1: 4, -22: 1}) * Fr(13, 19),
    '58a': Rad({1: 6, -22: 1}) * Fr(-5, 29),
    '58b': Rad({1: 6, -22: 1}) * Fr(3, 29),
}


# =========================================================================================== embedding and tau
def pq_coords(a, b, q):
    """Coordinates w.r.t. (1, P, Q, PQ) with P^2 = p > 0 (diagonal) and Q^2 = -qq < 0 (anti-diagonal)."""
    x, y, z, w = (Fr(t) for t in q)
    return (x, y, z, w, a, -b) if a > 0 else (x, z, y, -w, b, -a)


def matrix(a, b, q):
    """Image of q under the embedding (entries as Decimals)."""
    x, y, z, w, p, qq = pq_coords(a, b, q)
    sp, sq = Decimal(p).sqrt(), Decimal(qq).sqrt()
    dec = lambda f: Decimal(f.numerator) / Decimal(f.denominator)
    x, y, z, w = map(dec, (x, y, z, w))
    return (x + y * sp, (z + w * sp) * sq, -(z - w * sp) * sq, x - y * sp)


def tau_exact(a, b, mu, D):
    """tau = (y sqrt p + i s sqrt|D|/2) / C,  C = (w sqrt p - z) sqrt q,  s = sign C (so that Im tau > 0);
    rationalised with 1/C = (w sqrt p + z) sqrt q / ((w^2 p - z^2) q)."""
    x, y, z, w, p, qq = pq_coords(a, b, mu)
    dec = lambda f: Decimal(f.numerator) / Decimal(f.denominator)
    s = 1 if dec(w) * Decimal(p).sqrt() - dec(z) > 0 else -1
    num = Rad({p: y, D: Fr(s, 2)}) * Rad({p: w, 1: z}) * Rad({qq: 1})
    return num * Fr(1, (w * w * p - z * z) * qq)


def fixed_point(M):
    A, B, C, D_ = M
    disc = (D_ - A) ** 2 + 4 * B * C          # < 0 for an elliptic element
    re = (A - D_) / (2 * C)
    im = (-disc).sqrt() / (2 * abs(C))
    return re, im


def XY(M, tau):
    """Nelson's X, Y of sigma_tau^{-1} M sigma_tau, sigma_tau = [[y^1/2, x y^-1/2], [0, y^-1/2]]."""
    x, y = tau
    s1, s4 = y.sqrt(), 1 / y.sqrt()
    s2 = x * s4
    a, b, c, d = M
    A_, B_, C_, D_ = s4 * a - s2 * c, s4 * b - s2 * d, s1 * c, s1 * d
    P, Q, R, S = A_ * s1, A_ * s2 + B_ * s4, C_ * s1, C_ * s2 + D_ * s4
    return ((P + S) / 2, (Q - R) / 2), ((P - S) / 2, (Q + R) / 2)


# =========================================================================================== W_k (trapezoid rule)
def _stirling2(n):
    S = [[0] * (n + 2) for _ in range(n + 2)]
    S[0][0] = 1
    for i in range(1, n + 2):
        for k in range(1, i + 1):
            S[i][k] = k * S[i - 1][k] + S[i - 1][k - 1]
    return S


_S2 = _stirling2(20)
LICOEF = {s: [math.factorial(k) * _S2[s + 1][k + 1] for k in range(s + 1)] for s in range(15)}


def _choose_h_T(x, numax, digits):
    """Step and range for W_k(x) = 2^{1-k} int_0^oo [x cosh((k-1)t) Li_{-k-1}(q) - cosh(kt) Li_{-k}(q)] dt,
    q = exp(-x cosh t): the integrand is analytic in |Im t| < pi/2, so the trapezoid error is
    ~ exp(x (1 - cos a) - 2 pi a / h) relative to e^{-x}; we use a = 3/2."""
    a = 1.5
    ca = math.cos(a)
    xa = x * ca
    rho = math.exp(-xa / 2)
    lmbar = math.log(2 * math.gamma(numax) * 2.0 ** (2 * numax - 1)) - numax * math.log(xa) + \
        math.log(x * rho / (1 - rho) ** 2 + rho / (1 - rho))
    lmbar = min(lmbar, math.log(2 * (x + 1)) + 0.5 * math.log(2 * math.pi / xa) + numax ** 2 / (2 * xa) - xa + 30)
    log_eps = -digits * math.log(10) - x
    h = 2 * math.pi * a / (lmbar - log_eps + 2)
    T = 0.0
    while math.log(4 / 3 * (x + 1)) + numax * T - x * math.cosh(T) > log_eps - 10 or x * math.cosh(T) < 30:
        T += 0.02
    return Decimal('%.8g' % h), int(T / h) + 2


def W_values(x, nus, digits=WDIG):
    """{k: W_k(x)} for k in nus (x a Decimal > 0)."""
    numax = max(nus)
    hstep, J = _choose_h_T(float(x), numax, digits)
    acc = {k: Decimal(0) for k in nus}
    for jn in range(J + 1):
        t = jn * hstep
        et = t.exp()
        emt = 1 / et
        q = (-x * (et + emt) / 2).exp()
        r = q / (1 - q)
        rp = [Decimal(1)]
        for _ in range(numax + 2):
            rp.append(rp[-1] * r)
        Li = {s: sum(c * rp[k + 1] for k, c in enumerate(LICOEF[s])) for s in range(min(nus), numax + 2)}
        ep, em = [Decimal(1)], [Decimal(1)]
        for _ in range(numax):
            ep.append(ep[-1] * et)
            em.append(em[-1] * emt)
        ch = [(ep[m] + em[m]) / 2 for m in range(numax + 1)]
        wgt = Decimal(1) / 2 if jn == 0 else Decimal(1)
        for k in nus:
            acc[k] += wgt * (x * ch[k - 1] * Li[k + 1] - ch[k] * Li[k])
    return {k: acc[k] * hstep / Decimal(2) ** (k - 1) for k in nus}


# =========================================================================================== constants
def pi_dec():
    """pi by Machin's formula (working precision + 10 digits)."""
    from decimal import localcontext
    with localcontext() as ctx:
        ctx.prec = DIGITS + 10
        eps = Decimal(10) ** (-(DIGITS + 10))

        def atan_inv(n):
            s, term, k = Decimal(0), Decimal(1) / n, 0
            while term > eps:
                s += term / (2 * k + 1) if k % 2 == 0 else -term / (2 * k + 1)
                term /= n * n
                k += 1
            return s
        v = 16 * atan_inv(5) - 4 * atan_inv(239)
    return +v


def pi_omega2(D, PI):
    """pi Omega_D^2 = pi |eta((1+sqrt D)/2)|^4 = pi e^{-pi sqrt|D|/6} prod (1 - q^n)^4,  q = -e^{-pi sqrt|D|}."""
    s = Decimal(-D).sqrt()
    q = -(-PI * s).exp()
    prod, qn = Decimal(1), q
    while abs(qn) > Decimal(10) ** (-DIGITS - 5):
        prod *= (1 - qn)
        qn *= q
    return PI * (-PI * s / 6).exp() * prod ** 4


def omega2_gamma_float(D):
    """The paper's definition Omega_D = (2 pi |D|)^{-1/2} (prod Gamma(j/|D|)^{chi(j)})^{1/2}, double precision."""
    p = -D
    lg = sum((1 if pow(j, (p - 1) // 2, p) == 1 else -1) * math.lgamma(j / p) for j in range(1, p))
    return math.exp(lg) / (2 * math.pi * p)


# =========================================================================================== the CM computation
class CMPoint:
    def __init__(self, N):
        dd = DATA[N]
        self.N, self.a, self.b, self.D = N, dd['a'], dd['b'], dd['D']
        A = self.A = QuatAlg(dd['a'], dd['b'])
        mu = self.mu = tuple(Fr(t) for t in dd['mu'])
        if dd['gens'] is None:          # the unique maximal order containing Z<i, j, mu>
            orders = maximal_orders_containing(A, [(0, 1, 0, 0), (0, 0, 1, 0), mu], N)
            assert len(orders) == 1
            self.R = orders[0]
        else:
            self.R = ring_closure(A, dd['gens'])
        assert reduced_discriminant(A, self.R) == N and contains(self.R, mu)
        assert 2 * mu[0] == 1 and 1 - 4 * A.nrd(mu) == self.D
        self.p = -self.D
        one = (Fr(1), Fr(0), Fr(0), Fr(0))
        r_el = (2 * mu[0] - 1, 2 * mu[1], 2 * mu[2], 2 * mu[3])        # r = 2 mu - 1, r^2 = D
        # j0: a nonzero element orthogonal (trace form) to 1 and mu
        sol = None
        for c in itertools.product(range(-3, 4), repeat=3):
            v = (Fr(0), Fr(c[0]), Fr(c[1]), Fr(c[2]))
            if any(c) and A.pair(v, mu) == 0 and (sol is None or abs(A.nrd(v)) < abs(A.nrd(sol))):
                sol = v
        if sol is None:                 # general solution of trd(v)=0, trd(v conj(mu))=0
            a_, b_ = A.a, A.b
            u = (-a_ * mu[1], -b_ * mu[2], a_ * b_ * mu[3])         # pair(v,mu)/2 = v.u for v = (0, v1, v2, v3)
            e = [(Fr(0), -u[1], u[0], Fr(0)), (Fr(0), -u[2], Fr(0), u[0]), (Fr(0), Fr(0), -u[2], u[1])]
            sol = next(v for v in e if any(v))
        self.j0 = sol
        self.c = -A.nrd(sol)
        assert self.c > 0
        T = [list(one), list(r_el), list(sol), list(A.mul(r_el, sol))]
        self.Ti = inv(T)
        # numerics of the embedding
        Mmu = matrix(self.a, self.b, mu)
        self.tau = fixed_point(Mmu)
        Xr, _ = XY(matrix(self.a, self.b, r_el), self.tau)
        self.eps = 1 if Xr[1] > 0 else -1
        X0, self.Y0 = XY(matrix(self.a, self.b, sol), self.tau)
        assert abs(X0[0]) + abs(X0[1]) < Decimal(10) ** (-DIGITS + 10)
        assert abs(self.Y0[0] ** 2 + self.Y0[1] ** 2 - Decimal(self.c.numerator) / self.c.denominator) < Decimal(10) ** (-DIGITS + 10)
        self.tau_exact = tau_exact(self.a, self.b, mu, self.D)
        Rs = dual(A, self.R)
        self.divisors = [d for d in range(1, N + 1) if N % d == 0]
        self.split = {}
        for d in self.divisors:
            basis = lattice_basis(list(self.R) + [tuple((N // d) * x for x in v) for v in Rs])
            C = [[sum(v[k] * self.Ti[k][j] for k in range(4)) for j in range(4)] for v in basis]
            L = reduce(lcm, [x.denominator for row in C for x in row], 1)
            H = hnf_rows([[int(x * L) for x in row] for row in C])
            assert len(H) == 4 and H[1][0] == 0 and H[2][:2] == [0, 0] and H[3][:3] == [0, 0, 0]
            self.split[d] = (L, H)

    def grouped_sums(self, d, an, Mcut, nmax):
        """Exact G_n(d, K1) = sum a_f(d nrd beta) Z1^{2+n} Z2^n over beta in R^(d), L^2 N(kappa1) = K1, nrd beta > 0,
        Z1 = L kappa1 = X0 + X1 r, Z2 = L conj(kappa2) = X2 - X3 r (elements of Z[r], r^2 = D, as integer pairs)."""
        L, H = self.split[d]
        p, D = self.p, self.D
        cn, cd = self.c.numerator, self.c.denominator
        (h11, h12, h13, h14), (_, h22, h23, h24), (_, _, h33, h34), h44 = H[0], H[1], H[2], H[3][3]
        K1max = int(Mcut * L * L)

        def powers(P, Q, start, count):
            re, im, out = 1, 0, []
            for _ in range(start):
                re, im = re * P + im * Q * D, re * Q + im * P
            for _ in range(count):
                out.append((re, im))
                re, im = re * P + im * Q * D, re * Q + im * P
            return out
        byclass = defaultdict(lambda: defaultdict(lambda: [[0, 0] for _ in range(nmax + 1)]))
        l1max = math.isqrt(K1max) // h11 + 1
        for l1 in range(-l1max, l1max + 1):
            X0 = l1 * h11
            rem = K1max - X0 * X0
            if rem < 0:
                continue
            t = math.isqrt(rem // p) + 1
            for l2 in range((-t - l1 * h12) // h22 - 1, (t - l1 * h12) // h22 + 2):
                X1 = l1 * h12 + l2 * h22
                K1 = X0 * X0 + p * X1 * X1
                if K1 == 0 or K1 > K1max:
                    continue
                s2, s3 = l1 * h13 + l2 * h23, l1 * h14 + l2 * h24
                qq = s2 // h33
                s2, s3 = s2 - qq * h33, (s3 - qq * h34) % h44
                acc = byclass[(s2, s3)][K1]
                for n, (u, v) in enumerate(powers(X0, X1, 2, nmax + 1)):
                    acc[n][0] += u
                    acc[n][1] += v
        G = defaultdict(lambda: [[0, 0] for _ in range(nmax + 1)])
        den = cd * L * L
        for (s2, s3), dct in byclass.items():
            K2max = (cd * max(dct) - 1) // cn
            A2 = defaultdict(lambda: [[0, 0] for _ in range(nmax + 1)])
            m1max = math.isqrt(K2max) // h33 + 2
            for m1 in range(-m1max, m1max + 1):
                X2 = s2 + m1 * h33
                rem = K2max - X2 * X2
                if rem < 0:
                    continue
                base = s3 + m1 * h34
                t = math.isqrt(rem // p) + 1
                for m2 in range((-t - base) // h44 - 1, (t - base) // h44 + 2):
                    X3 = base + m2 * h44
                    K2 = X2 * X2 + p * X3 * X3
                    if K2 > K2max:
                        continue
                    acc = A2[K2]
                    for n, (u, v) in enumerate(powers(X2, -X3, 0, nmax + 1)):
                        acc[n][0] += u
                        acc[n][1] += v
            K2s = sorted(A2)
            for K1, J in dct.items():
                I = [[0, 0] for _ in range(nmax + 1)]
                for K2 in K2s:
                    num = cd * K1 - cn * K2
                    if num <= 0:
                        break
                    dn, rr = divmod(d * num, den)
                    assert rr == 0
                    av = an[dn]
                    if av:
                        row = A2[K2]
                        for n in range(nmax + 1):
                            I[n][0] += av * row[n][0]
                            I[n][1] += av * row[n][1]
                g = G[K1]
                for n in range(nmax + 1):
                    (jr, ji), (ir, ii) = J[n], I[n]
                    g[n][0] += jr * ir + D * ji * ii
                    g[n][1] += jr * ii + ji * ir
        return L, G

    def run(self, form, nmax, Mcut=MCUT, PI=None):
        PI = PI or pi_dec()
        N = self.N
        bad = [q for q in range(2, N + 1) if N % q == 0 and all(q % r for r in range(2, q))]
        an = an_table(CURVES[form], N * Mcut + 10, bad)
        sp = Decimal(self.p).sqrt()
        Sre, Sim = [Decimal(0)] * (nmax + 1), [Decimal(0)] * (nmax + 1)
        Wcache = {}
        for d in self.divisors:
            L, G = self.grouped_sums(d, an, Mcut, nmax)
            cf = Fr(1, d * an[d])
            for K1, g in G.items():
                m = Fr(K1, L * L)
                if m not in Wcache:
                    x = 4 * PI * (Decimal(m.numerator) / Decimal(m.denominator)).sqrt()
                    Wcache[m] = W_values(x, list(range(2, 3 + nmax)))
                Wv = Wcache[m]
                sK = Decimal(K1).sqrt()
                base = Decimal(cf.numerator) / Decimal(cf.denominator) / K1
                for n in range(nmax + 1):
                    gr, gi = g[n]
                    if gr or gi:
                        fac = base / (Decimal(L) ** n * sK ** n) * Wv[2 + n]
                        Sre[n] += gr * fac
                        Sim[n] += self.eps * gi * sp * fac          # A + B r  ->  A + B eps i sqrt|D|
        S = [complex_(Sre[n], Sim[n]) for n in range(nmax + 1)]
        c = {n: cdiv(cmul(S[n], cpow(S[0], n - 1)), cpow(S[1], n)) for n in range(2, nmax + 1)}
        Theta = cmul(cmul((-4 * PI, Decimal(0)), self.Y0), cdiv(S[1], S[0]))
        return c, Theta


# small complex helpers on (re, im) pairs of Decimals
def complex_(re, im):
    return (re, im)


def cmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def cdiv(u, v):
    n = v[0] * v[0] + v[1] * v[1]
    return ((u[0] * v[0] + u[1] * v[1]) / n, (u[1] * v[0] - u[0] * v[1]) / n)


def cpow(u, k):
    r = (Decimal(1), Decimal(0))
    for _ in range(k):
        r = cmul(r, u)
    return r


def cabs(u):
    return (u[0] * u[0] + u[1] * u[1]).sqrt()


def fmt(x, nd=45):
    return format(x, '.%dg' % nd)


def main(levels):
    PI = pi_dec()
    worst_c, worst_t = Decimal(0), Decimal(0)
    for N in levels:
        t0 = time.time()
        E = CMPoint(N)
        P = pi_omega2(E.D, PI)
        rel = abs(float(P) / (math.pi * omega2_gamma_float(E.D)) - 1)
        print('=' * 110)
        print('X_%d   B = (%d, %d),  mu = %s,  D = %d,  j0 = %s,  -nrd(j0) = %s' % (
            N, E.a, E.b, tuple(str(x) for x in E.mu), E.D, tuple(str(x) for x in E.j0), E.c))
        print('   tau_N (exact, from the data) = %s' % E.tau_exact)
        print('   tau_N in Table 5             = %s      identical: %s' % (TABLE5_TAU[N], E.tau_exact == TABLE5_TAU[N]))
        tr, ti = E.tau_exact.value()
        print('   tau_N = %s + %s i   (|exact - fixed point| = %.1e)' % (
            fmt(tr, 40), fmt(ti, 40), float(abs(tr - E.tau[0]) + abs(ti - E.tau[1]))))
        print('   pi Omega_D^2 = %s   (eta product; vs Gamma product in double precision: rel. diff %.1e)' % (fmt(P, 40), rel))
        for form in DATA[N]['forms']:
            nmax = max(TABLE4[form])
            c, Theta = E.run(form, nmax, PI=PI)
            print('   --- %s' % form)
            for n in range(2, nmax + 1):
                ex = TABLE4[form][n]
                exd = Decimal(ex.numerator) / Decimal(ex.denominator)
                diff = cabs((c[n][0] - exd, c[n][1]))
                worst_c = max(worst_c, diff)
                print('      c_%d = %-24s computed %s  |c_num - c_table| = %.1e' % (n, ex, fmt(c[n][0]), float(diff)))
            alg = TABLE5_THETA[form]
            ar, ai = alg.value()
            tab = cmul((P, Decimal(0)), (ar, ai))
            diff = cabs((Theta[0] - tab[0], Theta[1] - tab[1]))
            worst_t = max(worst_t, diff)
            cfmt = lambda u: '%s %s %s i' % (fmt(u[0]), '-' if u[1] < 0 else '+', fmt(abs(u[1])))
            print('      Theta_N = b_1/b_0 computed      = %s' % cfmt(Theta))
            print('      Table 5: pi Omega_D^2 * alpha_N = %s,   alpha_N = %s' % (cfmt(tab), alg))
            print('      |Theta_num - Theta_table| = %.1e' % float(diff))
        print('   (%.0f s)' % (time.time() - t0), flush=True)
    print('=' * 110)
    print('max |c_num - c_table| = %.1e     max |Theta_num - Theta_table| = %.1e' % (float(worst_c), float(worst_t)))


def list_orders(N):
    dd = DATA[N]
    A = QuatAlg(dd['a'], dd['b'])
    mu = tuple(Fr(t) for t in dd['mu'])
    orders = maximal_orders_containing(A, [(0, 1, 0, 0), (0, 0, 1, 0), mu], N)
    print('X_%d: %d maximal order(s) containing Z<i, j, mu>' % (N, len(orders)))
    PI = pi_dec()
    for O in orders:
        DATA[N] = dict(dd, gens=O)
        E = CMPoint(N)
        c, Theta = E.run(dd['forms'][0], 2, Mcut=36, PI=PI)
        P = pi_omega2(E.D, PI)
        print('   Z-basis', [tuple(str(x) for x in v) for v in O])
        print('      Theta/(pi Omega^2) = %.15f %+.15f i  (Table 5: %s);  c_2 = %.15f' % (
            float(Theta[0] / P), float(Theta[1] / P), TABLE5_THETA[dd['forms'][0]], float(c[2][0])))
    DATA[N] = dd


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == '--orders':
        for N in map(int, args[1:]):
            list_orders(N)
    else:
        main([int(x) for x in args] if args else [14, 15, 21, 33, 34, 46, 26, 38, 58])
