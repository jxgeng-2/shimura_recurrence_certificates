"""CM integration, numerical recognition and exact point lifting.
Loaded by Compute Phi and CM Images.sage; no computation starts on loading.
"""
# common/heegner_image_from_series.sage
"""CM divisor images from a disk expansion, with reusable Dirichlet reduction.

The side-data JSON is read without modifying it. Exact quaternion coordinates
are preferred to exported decimal matrices. Endpoint reduction computes an
integral modulo source periods; mu must send these periods into the target
elliptic-curve lattice. Use method='path' for actual path integrals/periods.

    reducer = DirichletReducer(side_data(21), prec_bits=1000)
    ans = heegner_image_from_series('w3', tauD, b=b, tau0=reducer.tau,
        AL=AL, mu=mu, E_label='21a1', prec_bits=1000, Nper=1000,
        reducer=reducer)

Tail estimates and identity detection are numerical diagnostics, not certified
error enclosures. No claim of algebraic recognition is made by this module.
"""

from sage.all import (QQ, ZZ, ComplexField, RealField, EllipticCurve,
                      matrix, vector, identity_matrix, log)
import json
import os
import re


def _as_CC_matrix(g, CC):
    try:
        return matrix(CC, 2, 2, [CC(g[0,0]), CC(g[0,1]), CC(g[1,0]), CC(g[1,1])])
    except (TypeError, IndexError):
        return matrix(CC, g)


def _mobius(g, z):
    return (g[0,0]*z + g[0,1])/(g[1,0]*z + g[1,1])


def _inv_cayley(w, p, pbar):
    return (pbar*w-p)/(w-1)


def _cayley(z, p, pbar):
    return (z-p)/(z-pbar)


def _apply_disk(g, w, p, pbar):
    return _cayley(_mobius(g, _inv_cayley(w, p, pbar)), p, pbar)


def _heegner_quaternion_matrix(coords, a, b, RF):
    x, y, z, t = [QQ(c) for c in coords]
    a, b = QQ(a), QQ(b)
    if a <= 0 and b > 0:
        a, b, y, z, t = b, a, z, y, -t
    if a <= 0:
        raise ValueError('The quaternion algebra needs a positive real Hilbert parameter.')
    sa, sb = RF(a).sqrt(), RF(abs(b)).sqrt()
    sign = 1 if b > 0 else -1
    return matrix(RF, [[RF(x)+RF(y)*sa, (RF(z)+RF(t)*sa)*sb],
                      [sign*(RF(z)-RF(t)*sa)*sb, RF(x)-RF(y)*sa]])


def _heegner_fixed_point(g, CC):
    a, b, c, d = [CC(t) for t in g.list()]
    disc = (a-d)**2 + 4*b*c
    roots = [(a-d+disc.sqrt())/(2*c), (a-d-disc.sqrt())/(2*c)]
    return next(z for z in roots if z.imag() > 0)


def _heegner_read_side_data(source):
    if isinstance(source, dict):
        return dict(source)
    with open(os.fspath(source), 'r', encoding='utf-8-sig') as stream:
        text = stream.read()
    # Some Magma exports contain presentation-only line continuations.
    text = re.sub(r'\\\r?\n[ \t]*', '', text)
    text = re.sub(r'\\([\[\]_])', r'\1', text)
    return json.loads(text)


class DirichletReducer:
    """Repeated side pairing using hyperbolic distance from the disk center.

    A Dirichlet domain centered at tau is characterized by |g.w| >= |w|
    for all its side pairings. Each accepted side strictly lowers |w|.
    The returned gamma acts on the upper half-plane, not on disk coordinates.
    """
    INPUT_API_VERSION = 2

    def __init__(self, side_data=None, prec_bits=300, tau0=None,
                 S_matrices=None, max_steps=10000):
        self.prec_bits = int(prec_bits)
        self.CC = ComplexField(self.prec_bits)
        self.RF = RealField(self.prec_bits)
        self.max_steps = int(max_steps)
        self.data = None if side_data is None else _heegner_read_side_data(side_data)
        CC, RF = self.CC, self.RF
        if self.data is not None:
            data = self.data
            if all(k in data for k in ('cm_element', 'quaternion_a', 'quaternion_b')):
                a, b = QQ(data['quaternion_a']), QQ(data['quaternion_b'])
                # Expose exact geometry for the companion export's make_reducer.
                self.a, self.b = a, b
                self.cm_element = tuple(QQ(c) for c in data['cm_element'])
                self.side_quaternions = [tuple(QQ(c) for c in row)
                                         for row in data['side_generators']]
                center = _heegner_fixed_point(
                    _heegner_quaternion_matrix(data['cm_element'], a, b, RF), CC)
                generators = []
                for coordinates in data['side_generators']:
                    x, y, z, t = [QQ(c) for c in coordinates]
                    norm = x*x-a*y*y-b*z*z+a*b*t*t
                    if norm != 1:
                        raise ValueError('A side pairing has exact norm %s, expected 1.' % norm)
                    generators.append(_heegner_quaternion_matrix(coordinates, a, b, RF))
                self.exact_coordinates = True
            else:
                center = CC(RF(data['tau_re']), RF(data['tau_im']))
                generators = [matrix(RF, [[RF(c) for c in row] for row in g])
                              for g in data['side_matrices']]
                self.exact_coordinates = False
            if tau0 is not None:
                supplied = CC(tau0)
                tolerance = RF(2)**(-min(80, self.prec_bits//3))
                if abs(supplied-center) > tolerance*max(1, abs(center)):
                    raise ValueError('tau0 differs from the center in the side-data JSON.')
                center = supplied
        else:
            if tau0 is None or S_matrices is None:
                raise ValueError('Supply side_data, or tau0 together with S_matrices.')
            center = CC(tau0)
            generators = list(S_matrices)
            self.exact_coordinates = False
        self.tau = CC(center)
        self.tau0 = self.tau
        self.pbar = self.tau.conjugate()
        if self.tau.imag() <= 0:
            raise ValueError('The expansion center must be in the upper half-plane.')
        self.tol = RF(2)**(-max(20, self.prec_bits-40))
        self.generators = [_as_CC_matrix(g, CC) for g in generators]
        # Include inverse pairings, and deduplicate only for speed. The actual
        # matrices retain full precision; rounded keys never enter arithmetic.
        C = matrix(CC, [[1, -self.tau], [1, -self.pbar]])
        Cinv = C.inverse()
        self.actions = []
        seen = set()
        for g in self.generators + [g.inverse() for g in self.generators]:
            if g.det().real() <= 0:
                raise ValueError('Reduction matrices must preserve the upper half-plane.')
            flat = list(g.list())
            scale = next(t for t in flat if abs(t) > self.tol)
            key = tuple((round(float((t/scale).real()), 12),
                         round(float((t/scale).imag()), 12)) for t in flat)
            if key in seen:
                continue
            seen.add(key)
            disk = C*g*Cinv
            self.actions.append((g, disk))

    def reduce(self, w, return_metadata=False):
        CC = self.CC
        w_original = CC(w)
        current = w_original
        if abs(current) >= 1:
            raise ValueError('The point to reduce must be in the open unit disk.')
        gamma = identity_matrix(CC, 2)
        chosen = []
        for count in range(self.max_steps):
            radius = abs(current)
            best = None
            best_w = current
            for index, (g, disk) in enumerate(self.actions):
                image = _mobius(disk, current)
                image_radius = abs(image)
                if image_radius < radius - self.tol:
                    best, best_w, radius = index, image, image_radius
            if best is None:
                break
            gamma = self.actions[best][0]*gamma
            current = best_w
            chosen.append(best)
        else:
            raise ArithmeticError('Side reduction exceeded max_steps=%s.' % self.max_steps)
        metadata = {'steps': len(chosen), 'side_word': chosen,
                    'initial_radius': abs(w_original), 'reduced_radius': abs(current),
                    'num_sides': len(self.actions)}
        if return_metadata:
            return gamma, current, metadata
        return gamma, current

    def reduce_upper(self, z, return_metadata=False):
        return self.reduce(_cayley(self.CC(z), self.tau, self.pbar), return_metadata)


def reduce_point_disc(w, reducer, prec_bits=None, return_metadata=False):
    """Notebook-compatible wrapper; reuse the reducer's precision and center."""
    if prec_bits is not None and int(prec_bits) != reducer.prec_bits:
        raise ValueError('Construct the reducer at the requested precision.')
    return reducer.reduce(w, return_metadata=return_metadata)


def reduction_words(Smats, depth=3):
    """Legacy word enumeration, retained for callers without side data."""
    def key(g):
        flat = list(g.list())
        scale = max(flat, key=lambda t: abs(t))
        return tuple(round(float((x/scale).real()), 10) for x in flat)
    out, seen, level = [], set(), list(Smats)
    for unused in range(int(depth)):
        nxt = []
        for g in level:
            k = key(g)
            if k not in seen:
                seen.add(k)
                out.append(g)
                nxt.extend(s*g for s in Smats)
        level = nxt
    return out


def prepare_reduction(tau0=None, S_matrices=None, prec_bits=300,
                      reduction_depth=3, side_data=None):
    if side_data is not None:
        return DirichletReducer(side_data, prec_bits=prec_bits, tau0=tau0)
    CC = ComplexField(prec_bits)
    matrices = [_as_CC_matrix(g, CC) for g in S_matrices]
    words = reduction_words(matrices, reduction_depth)
    return DirichletReducer(prec_bits=prec_bits, tau0=tau0, S_matrices=words)


def reduce_disk(w, words, p, pbar, CC, tol=None):
    """Legacy signature, accepting either words or a prepared reducer."""
    reducer = words if isinstance(words, DirichletReducer) else DirichletReducer(
        prec_bits=CC.precision(), tau0=p, S_matrices=words)
    return reducer.reduce(w)


def _heegner_resolve_reducer(tau0, S_matrices, prec_bits, reduction_depth,
                            reducer=None, side_data=None):
    if reducer is None:
        reducer = prepare_reduction(tau0, S_matrices, prec_bits,
                                    reduction_depth, side_data)
    if reducer.prec_bits != int(prec_bits):
        raise ValueError('Reducer precision and integration precision differ.')
    RF = reducer.RF
    if abs(reducer.CC(tau0)-reducer.tau) > RF(2)**(-min(80, int(prec_bits)//3)):
        raise ValueError('Reducer and power series use different centers.')
    return reducer


def antiderivative_evaluator(b, tau0, Nper=None, prec_bits=300):
    """H(u)=(tau0-conj(tau0))*sum b[n]*u^(n+1)/(n+1), by Horner."""
    CC, RF = ComplexField(prec_bits), RealField(prec_bits)
    p = CC(tau0)
    Nper = len(b) if Nper is None else int(Nper)
    if Nper < 1 or Nper > len(b):
        raise ValueError('Nper must be between 1 and the available coefficient count.')
    h = [CC(0)] + [(p-p.conjugate())*CC(b[n])/CC(n+1) for n in range(Nper)]
    def H(u):
        u, value = CC(u), CC(0)
        for coefficient in reversed(h[1:]):
            value = (value+coefficient)*u
        return value
    def tail(u):
        # Sum absolute contributions in the final two blocks. Extrapolation
        # uses their observed decay; this is explicitly an empirical estimate.
        radius = RF(abs(CC(u)))
        if radius == 0:
            return RF(0)
        width = min(32, max(1, Nper//4))
        end = sum(RF(abs(h[n]))*radius**n
                  for n in range(max(1,Nper-width+1), Nper+1))
        prior = sum(RF(abs(h[n]))*radius**n
                    for n in range(max(1,Nper-2*width+1), max(1,Nper-width+1)))
        ratio = end/prior if prior else RF(1)
        return end/max(RF('0.01'), 1-ratio) if ratio < 1 else end*100
    H.tail_estimate = tail
    H.Nper = Nper
    return H


def _heegner_stats(reducer, reduction_depth, max_radius):
    return {'segments': 0, 'bad_segments': 0, 'max_u': 0.0,
            'num_words': len(reducer.actions), 'num_sides': len(reducer.actions),
            'reduction_depth': int(reduction_depth), 'max_radius': float(max_radius),
            'tail_abs_estimate': reducer.RF(0), 'tail_is_rigorous': False}


def _heegner_add_segment(u0, u1, H, stats, max_radius):
    radius = max(abs(u0), abs(u1))
    stats['segments'] += 1
    stats['max_u'] = max(stats['max_u'], float(radius))
    stats['bad_segments'] += int(radius > max_radius)
    stats['tail_abs_estimate'] += H.tail_estimate(u0)+H.tail_estimate(u1)
    return H(u1)-H(u0)


def integrate_between_points_from_series(z0, z1, b, tau0, S_matrices=None,
                                         Nper=None, prec_bits=300,
                                         reduction_depth=3,
                                         initial_subdivisions=20,
                                         max_radius=0.93, max_subdivide=20,
                                         verbose=True, reducer=None, side_data=None):
    """Actual path integral, with side reduction of each common subsegment."""
    reducer = _heegner_resolve_reducer(tau0, S_matrices, prec_bits,
                                      reduction_depth, reducer, side_data)
    CC, p, pbar = reducer.CC, reducer.tau, reducer.pbar
    H = antiderivative_evaluator(b, p, Nper, prec_bits)
    w0, w1 = _cayley(CC(z0),p,pbar), _cayley(CC(z1),p,pbar)
    stats = _heegner_stats(reducer, reduction_depth, max_radius)
    stats['mode'] = 'path_side_reduction'
    def segment(wa, wb, depth):
        gamma, unused = reducer.reduce((wa+wb)/2)
        u0, u1 = _apply_disk(gamma,wa,p,pbar), _apply_disk(gamma,wb,p,pbar)
        if max(abs(u0),abs(u1)) <= max_radius or depth <= 0:
            return _heegner_add_segment(u0,u1,H,stats,max_radius)
        midpoint = (wa+wb)/2
        return segment(wa,midpoint,depth-1)+segment(midpoint,wb,depth-1)
    K = max(1,int(initial_subdivisions))
    value = sum((segment(w0+(w1-w0)*CC(j)/K, w0+(w1-w0)*CC(j+1)/K,
                         int(max_subdivide)) for j in range(K)), CC(0))
    if verbose:
        print('path: %s segments, max |u| %.6f, empirical tail %s' %
              (stats['segments'], stats['max_u'], stats['tail_abs_estimate'].n(20)))
    return value, w0, w1, stats


def integrate_between_points_by_endpoint_reduction(z0, z1, b, tau0, S_matrices=None,
                                                   Nper=None, prec_bits=300,
                                                   reduction_depth=3,
                                                   max_radius=0.93, verbose=True,
                                                   reducer=None, side_data=None):
    """Integral modulo source periods via independently reduced endpoints."""
    reducer = _heegner_resolve_reducer(tau0,S_matrices,prec_bits,reduction_depth,
                                      reducer,side_data)
    CC, p, pbar = reducer.CC,reducer.tau,reducer.pbar
    H = antiderivative_evaluator(b,p,Nper,prec_bits)
    w0, w1 = _cayley(CC(z0),p,pbar), _cayley(CC(z1),p,pbar)
    g0,u0,m0 = reducer.reduce(w0,True)
    g1,u1,m1 = reducer.reduce(w1,True)
    stats = _heegner_stats(reducer,reduction_depth,max_radius)
    stats.update({'mode':'endpoint_reduction', 'gamma0':g0,'gamma1':g1,
                  'w0_raw_abs':float(abs(w0)),'w1_raw_abs':float(abs(w1)),
                  'u0_abs':float(abs(u0)),'u1_abs':float(abs(u1)),
                  'reduction0':m0,'reduction1':m1})
    value = _heegner_add_segment(u0,u1,H,stats,max_radius)
    if verbose:
        print('endpoint reduction: |u0| %.6f, |u1| %.6f, %s sides, empirical tail %s' %
              (abs(u0),abs(u1),len(reducer.actions),stats['tail_abs_estimate'].n(20)))
    return value,u0,u1,stats


def integrate_basepoint_to_W_basepoint(W,b,tau0,S_matrices=None,Nper=None,
                                       prec_bits=300,reduction_depth=3,K_sub=200,
                                       verbose=True,reducer=None,side_data=None):
    CC = ComplexField(prec_bits)
    return integrate_between_points_by_endpoint_reduction(
        tau0,_mobius(_as_CC_matrix(W,CC),CC(tau0)),b,tau0,S_matrices,
        Nper,prec_bits,reduction_depth,verbose=verbose,reducer=reducer,side_data=side_data)


def _resolve_W(W, AL=None):
    if isinstance(W,str):
        AL = globals().get('AL') if AL is None else AL
        if AL is None:
            raise ValueError('Supply AL when W is a string.')
        return AL[W]
    return W


def heegner_image_from_series(W,tau_prime,b=None,tau0=None,S_matrices=None,
                              AL=None,mu=None,E_label='15a1',Nper=None,
                              prec_bits=300,reduction_depth=3,
                              initial_subdivisions=20,max_radius=0.93,
                              max_subdivide=20,basepoint_fast=True,
                              method='endpoint',verbose=True,reducer=None,
                              side_data=None,identity_min_bits=80):
    """Return the image of [W*tau_prime]-[tau_prime] and numerical diagnostics.

    ``side_data`` can be a JSON path or dictionary. Reuse a ``reducer`` for
    batch calculations. No word enumeration is used when side data is given.
    Identity classification is numerical and can be disabled with None.
    """
    b = globals().get('b_n_seq_pristine') if b is None else b
    tau0 = globals().get('tau') if tau0 is None else tau0
    if tau0 is None and reducer is not None:
        tau0 = reducer.tau
    S_matrices = globals().get('S_matrices') if S_matrices is None else S_matrices
    mu = globals().get('mu') if mu is None else mu
    if b is None or tau0 is None or mu is None:
        raise ValueError('Supply b, tau0, and the differential scale mu.')
    reducer = _heegner_resolve_reducer(tau0,S_matrices,prec_bits,reduction_depth,
                                      reducer,side_data)
    CC,RF = reducer.CC,reducer.RF
    tp = CC(tau_prime)
    Wtp = _mobius(_as_CC_matrix(_resolve_W(W,AL),CC),tp)
    common = dict(Nper=Nper,prec_bits=prec_bits,reduction_depth=reduction_depth,
                  max_radius=max_radius,verbose=verbose,reducer=reducer)
    if method in ('endpoint','endpoint_reduction','fast','basepoint','old'):
        integral,w_start,w_end,stats = integrate_between_points_by_endpoint_reduction(
            tp,Wtp,b,tau0,S_matrices,**common)
    elif method in ('path','generic'):
        integral,w_start,w_end,stats = integrate_between_points_from_series(
            tp,Wtp,b,tau0,S_matrices,initial_subdivisions=initial_subdivisions,
            max_subdivide=max_subdivide,**common)
    else:
        raise ValueError("Unknown method; choose 'endpoint' or 'path'.")
    E = EllipticCurve(E_label)
    lattice = E.period_lattice()
    periods = [CC(t) for t in lattice.basis(prec=int(prec_bits))]
    zE = CC(mu)*integral
    period_matrix = matrix(RF,[[t.real() for t in periods],[t.imag() for t in periods]])
    coefficients = period_matrix.solve_right(vector(RF,[zE.real(),zE.imag()]))
    nearest = [ZZ(t.round()) for t in coefficients]
    z_reduced = zE-sum(nearest[j]*periods[j] for j in range(2))
    distance = max(abs(coefficients[j]-nearest[j]) for j in range(2))
    identity_bits = int(prec_bits) if distance == 0 else max(0,int((-log(distance,2)).floor()))
    identity = identity_min_bits is not None and identity_bits >= int(identity_min_bits)
    P = E.change_ring(CC)(0) if identity else lattice.elliptic_exponential(z_reduced)
    error = abs(CC(mu))*stats['tail_abs_estimate']
    stats['zE_tail_abs_estimate'] = error
    stats['tail_bits_estimate'] = int(prec_bits) if error == 0 else max(0,int((-log(error,2)).floor()))
    stats['period_distance'] = distance
    stats['identity_bits'] = identity_bits
    if verbose:
        print('zE =',zE.n(30))
        print('P =', 'O (numerical identity)' if identity else P)
    return {'tau_prime':tp,'W_tau_prime':Wtp,'I':integral,'zE':zE,
            'zE_reduced':z_reduced,'P':P,'w_start':w_start,'w_end':w_end,
            'is_identity':identity,'identity_bits':identity_bits,'stats':stats}


# common/cm_image_table_package.sage
"""Numerical algebraic recognition and exact lifting into class fields."""

from sage.all import (
    QQ, ZZ, ComplexField, EllipticCurve, Integer, NumberField, PolynomialRing,
    QuadraticField, RealField, latex, load, log, matrix, save,
)
from sage.arith.misc import algdep

import os
import re
import warnings


CM_IMAGE_TABLE_PACKAGE_VERSION = 1


def _cm_exact_rational(value):
    if isinstance(value, str):
        return QQ(value)
    try:
        parent = value.parent()
        if parent.is_exact():
            return QQ(value)
    except (AttributeError, TypeError, ValueError):
        pass
    if isinstance(value, int):
        return QQ(value)
    raise TypeError(
        "Quaternion coordinates and Hilbert-symbol parameters must be exact; "
        "use QQ values or strings instead of floating-point numbers."
    )


def _cm_fixed_point_upper_half_plane(M, CC):
    """Return the fixed point of M lying in the upper half-plane."""
    aa, bb, cc, dd = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    discriminant = (dd - aa) ** 2 + 4 * bb * cc
    imaginary_unit = CC(-1).sqrt()

    if cc == 0:
        raise ValueError("The supplied elliptic matrix has c=0.")
    root1 = (
        (aa - dd) + imaginary_unit * CC(-discriminant).sqrt()
    ) / (2 * cc)
    root2 = (
        (aa - dd) - imaginary_unit * CC(-discriminant).sqrt()
    ) / (2 * cc)
    roots = [z for z in (root1, root2) if z.imag() > 0]
    if not roots:
        raise ValueError("The supplied matrix has no fixed point in H.")
    return roots[0]


def _cm_sanitize_name(value):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))
    return text.strip("_.") or "curve"


def _cm_is_infinity(P):
    try:
        return bool(P.is_zero())
    except (AttributeError, TypeError):
        pass
    try:
        return bool(P[2] == 0)
    except (IndexError, TypeError):
        return False


def _cm_numerical_identity_from_coordinates(x, y, min_bits=80):
    r"""
    Recognize a numerical approximation to the identity of an elliptic curve.

    ``elliptic_exponential`` can return enormous affine coordinates instead of
    the exact projective point at infinity when its argument is merely
    *numerically* equal to a period.  At the identity, the standard local
    parameter ``t = -x/y`` satisfies

        x = t^-2 + O(t^-1),    y = -t^-3 + O(t^-2).

    We therefore require all three compatible estimates
    ``log2(|x|)/2``, ``log2(|y|)/3``, and ``-log2(|x/y|)`` to exceed
    ``min_bits``.  The returned bit count is a numerical diagnostic, not a
    proof that an arbitrary finite point with enormous height is the identity.
    """
    if x is None or y is None:
        return False, None

    try:
        precision = max(
            80,
            int(x.parent().precision()),
            int(y.parent().precision()),
        )
    except (AttributeError, TypeError, ValueError):
        precision = 80

    try:
        RRp = RealField(precision)
        abs_x = RRp(abs(x))
        abs_y = RRp(abs(y))
        if abs_x <= 1 or abs_y <= 1:
            return False, None

        abs_t = RRp(abs(x / y))
        if abs_t == 0:
            return True, Integer(precision)
        if abs_t >= 1:
            return False, None

        identity_bits = Integer(min(
            log(abs_x, 2) / 2,
            log(abs_y, 2) / 3,
            -log(abs_t, 2),
        ).floor())
        return identity_bits >= Integer(min_bits), identity_bits
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
        return False, None


def _cm_tail_usable_bits(stats, Nper, prec_bits):
    arithmetic_bits = max(Integer(0), Integer(prec_bits) - 50)
    if not isinstance(stats, dict):
        return arithmetic_bits
    radius = stats.get("max_u", stats.get("max_radius_seen"))
    if radius is None:
        return arithmetic_bits
    RRp = RealField(max(80, int(prec_bits)))
    radius = RRp(radius)
    if radius <= 0:
        return arithmetic_bits
    if radius >= 1:
        return Integer(0)
    tail_bits = Integer(((-Integer(Nper) * log(radius, 2)).floor()))
    tail_bits = max(Integer(0), tail_bits - 70)
    return min(arithmetic_bits, tail_bits)


def _cm_primitive_positive(R, polynomial):
    f = R(polynomial)
    if f == 0:
        return f
    content = ZZ(f.content())
    if content:
        f //= content
    if f.leading_coefficient() < 0:
        f = -f
    return f


def _cm_relative_residual_bits(f, value):
    RRp = RealField(value.parent().precision())
    absolute_value = RRp(abs(value))
    scale = RRp(0)
    for exponent, coefficient in enumerate(f.list()):
        scale += (
            abs(ZZ(coefficient))
            * max(RRp(1), absolute_value) ** exponent
        )
    residual = RRp(abs(f(value)))
    if residual == 0:
        return Integer(value.parent().precision())
    if scale == 0:
        return Integer(0)
    ratio = residual / scale
    if ratio >= 1:
        return Integer(0)
    return Integer((-log(ratio, 2)).floor())


def _cm_cluster(values, known_bits):
    if not values:
        return []
    CCp = values[0].parent()
    RRp = RealField(CCp.precision())
    comparison_bits = max(20, min(Integer(known_bits) // 3, 160))
    relative_tolerance = RRp(2) ** (-comparison_bits)
    representatives = []
    for value in values:
        if all(
            abs(value - other)
            > relative_tolerance * max(RRp(1), abs(value), abs(other))
            for other in representatives
        ):
            representatives.append(value)
    return representatives


def _cm_algdep_precisions(known_bits, cap, ladder):
    limit = min(Integer(known_bits), Integer(cap))
    values = [Integer(bits) for bits in ladder if Integer(bits) <= limit]
    if limit >= 60 and limit not in values:
        values.append(limit)
    return sorted(set(values), reverse=True)


def _cm_coefficient_height(f):
    return max([abs(ZZ(coefficient)) for coefficient in f.list()] or [ZZ(0)])


def _cm_stability_windows(known_bits, cap, ladder, radius):
    """
    Return consecutive precision windows, from highest to lowest.

    With the default radius 1, a center p is tested at p-1, p, and p+1.
    The top window is shifted down when necessary so that it stays below the
    available/capped precision.
    """
    radius = Integer(radius)
    if radius < 1:
        raise ValueError("algdep_stability_radius must be at least 1.")
    limit = min(Integer(known_bits), Integer(cap))
    centers = {
        Integer(bits) for bits in ladder
        if Integer(bits) - radius >= 60
        and Integer(bits) + radius <= limit
    }
    if limit - radius >= 60:
        centers.add(limit - radius)
    return [
        tuple(range(
            int(center - radius),
            int(center + radius) + 1,
        ))
        for center in sorted(centers, reverse=True)
    ]


def _cm_algdep_candidates_at_precision(
    selected_x,
    use_bits,
    R,
    max_degree,
    algdep_height_bound,
    algdep_min_residual_bits,
    algdep_residual_safety_bits,
):
    """Return every irreducible candidate surviving the full-value residual."""
    working_value = ComplexField(use_bits)(selected_x)
    candidates = {}
    errors = []
    for degree in range(1, int(max_degree) + 1):
        try:
            candidate = algdep(
                working_value,
                degree,
                use_bits=use_bits,
                height_bound=algdep_height_bound,
            )
        except (ValueError, ArithmeticError) as error:
            errors.append(
                "%s bits, degree %s: %s" % (use_bits, degree, error)
            )
            continue
        if candidate is None:
            continue

        candidate = _cm_primitive_positive(R, candidate)
        if candidate == 0 or not candidate.is_irreducible():
            continue
        residual_bits = _cm_relative_residual_bits(candidate, selected_x)
        required_bits = max(
            Integer(algdep_min_residual_bits),
            Integer(use_bits) - Integer(algdep_residual_safety_bits),
        )
        if residual_bits < required_bits:
            continue

        key = tuple(ZZ(coefficient) for coefficient in candidate.list())
        old = candidates.get(key)
        if old is None:
            candidates[key] = {
                "polynomial": candidate,
                "residual_bits": residual_bits,
                "required_bits": required_bits,
                "request_degrees": {Integer(degree)},
            }
        else:
            old["request_degrees"].add(Integer(degree))
            old["residual_bits"] = max(old["residual_bits"], residual_bits)
    return candidates, errors


def _cm_default_hilbert_class_field(D, h):
    """
    Construct H_D as an absolute number field for a fundamental D.

    For nonfundamental discriminants a ring-class-field constructor must be
    supplied with ``class_field_function``.
    """
    D = Integer(D)
    h = Integer(h)
    quadratic_field = QuadraticField(D, "sqrtD")
    if Integer(quadratic_field.discriminant()) != D:
        raise ValueError(
            "D=%s is not a fundamental discriminant; supply "
            "class_field_function(D,h) for its ring class field." % D
        )
    computed_h = Integer(quadratic_field.class_number())
    if computed_h != h:
        raise ValueError(
            "Magma supplied h=%s for D=%s, but Sage computes h=%s."
            % (h, D, computed_h)
        )
    relative_field = quadratic_field.hilbert_class_field("eta")
    absolute_field = relative_field.absolute_field("theta")
    expected_degree = Integer(2) * h
    if Integer(absolute_field.absolute_degree()) != expected_degree:
        raise ArithmeticError(
            "The constructed H_%s has absolute degree %s, expected %s."
            % (D, absolute_field.absolute_degree(), expected_degree)
        )
    return absolute_field


def _cm_candidate_in_class_field(polynomial, class_field):
    """Test exactly whether Q[x]/(polynomial) embeds in class_field."""
    if polynomial.degree() == 1:
        return True, Integer(1)
    candidate_field = NumberField(
        polynomial.change_ring(QQ), "candidate_root"
    )
    embeddings = candidate_field.embeddings(class_field)
    return bool(embeddings), Integer(len(embeddings))


def _cm_recognize_x(
    selected_item,
    orbit_items,
    D,
    h,
    R,
    algdep_use_bits_cap,
    algdep_use_bits_ladder,
    algdep_height_bound,
    algdep_min_residual_bits,
    algdep_stability_radius,
    algdep_residual_safety_bits,
    max_degree,
    class_field,
    class_field_error,
    require_class_field_containment,
    require_full_orbit,
    algdep_holdout_bits=30,
    algdep_min_stable_windows=2,
    algdep_window_gap_bits=20,
):
    # Consecutive-bit stability alone accepts many continued-fraction
    # approximants to transcendental numbers.  Reserve independent digits
    # of the original input, and require stability at separated windows.
    algdep_holdout_bits = Integer(algdep_holdout_bits)
    algdep_min_stable_windows = Integer(algdep_min_stable_windows)
    algdep_window_gap_bits = Integer(algdep_window_gap_bits)
    if algdep_holdout_bits < 20:
        raise ValueError("algdep_holdout_bits must be at least 20.")
    if algdep_min_stable_windows < 2:
        raise ValueError("algdep_min_stable_windows must be at least 2.")
    if algdep_window_gap_bits < 20:
        raise ValueError("algdep_window_gap_bits must be at least 20.")
    identity_count = sum(
        bool(item.get("is_infinity", False)) for item in orbit_items
    )
    if identity_count:
        if identity_count != len(orbit_items):
            return {
                "status": "inconsistent_identity_orbit",
                "polynomial": None,
                "known_bits": min(
                    Integer(item["usable_bits"]) for item in orbit_items
                ),
                "message": (
                    "%s of %s supplied CM representatives look like O. "
                    "A Galois orbit cannot mix O with finite points."
                    % (identity_count, len(orbit_items))
                ),
            }
        return {
            "status": "identity",
            "polynomial": None,
            "known_bits": min(
                Integer(item["usable_bits"]) for item in orbit_items
            ),
            "identity_bits": selected_item.get("identity_bits"),
            "identity_method": selected_item.get(
                "identity_method", "exact_or_cached"
            ),
            "message": (
                "The image is the elliptic-curve identity O; "
                "there is no finite x-coordinate to recognize."
            ),
        }

    selected_x = selected_item["x"]
    known_bits = min(
        Integer(item["usable_bits"]) for item in orbit_items
    )
    stability_windows = _cm_stability_windows(
        known_bits - algdep_holdout_bits,
        algdep_use_bits_cap,
        algdep_use_bits_ladder,
        algdep_stability_radius,
    )
    if not stability_windows:
        return {
            "status": "insufficient_tail_precision",
            "polynomial": None,
            "known_bits": known_bits,
            "message": (
                "There are not enough usable bits for a consecutive "
                "algdep stability window."
            ),
        }

    finite_x_values = [item["x"] for item in orbit_items]
    full_numeric_orbit = (
        finite_x_values + [value.conjugate() for value in finite_x_values]
    )

    candidates_by_precision = {}
    errors = []
    stable_candidates = {}
    rejected_field = []
    rejected_orbit = []
    rejected_holdout = []
    containment_cache = {}

    for window in stability_windows:
        for use_bits in window:
            if use_bits in candidates_by_precision:
                continue
            candidates, new_errors = _cm_algdep_candidates_at_precision(
                selected_x,
                Integer(use_bits),
                R,
                max_degree,
                algdep_height_bound,
                algdep_min_residual_bits,
                algdep_residual_safety_bits,
            )
            candidates_by_precision[use_bits] = candidates
            errors.extend(new_errors)

        common_keys = set(candidates_by_precision[window[0]])
        for use_bits in window[1:]:
            common_keys.intersection_update(
                candidates_by_precision[use_bits]
            )

        accepted_in_window = []
        for key in common_keys:
            common_request_degrees = set(
                candidates_by_precision[window[0]][key]["request_degrees"]
            )
            for use_bits in window[1:]:
                common_request_degrees.intersection_update(
                    candidates_by_precision[use_bits][key]["request_degrees"]
                )
            if not common_request_degrees:
                continue
            polynomial = candidates_by_precision[window[0]][key][
                "polynomial"
            ]
            stable_entry = stable_candidates.setdefault(key, {
                "polynomial": polynomial,
                "windows": [],
                "accepted_windows": [],
            })
            stable_entry["windows"].append(window)

            holdout_required_bits = max(
                Integer(algdep_min_residual_bits),
                Integer(max(window)) + algdep_holdout_bits,
            )
            selected_residual_bits = min(
                candidates_by_precision[use_bits][key]["residual_bits"]
                for use_bits in window
            )
            if selected_residual_bits < holdout_required_bits:
                rejected_holdout.append((
                    polynomial, selected_residual_bits, holdout_required_bits
                ))
                continue

            if require_class_field_containment:
                if class_field_error is not None or class_field is None:
                    rejected_field.append((
                        polynomial,
                        "class field unavailable: %s" % class_field_error,
                    ))
                    continue
                if key not in containment_cache:
                    try:
                        containment_cache[key] = (
                            _cm_candidate_in_class_field(
                                polynomial, class_field
                            )
                        )
                    except Exception as error:
                        containment_cache[key] = (False, Integer(0))
                        rejected_field.append((
                            polynomial,
                            "embedding computation failed: %s" % error,
                        ))
                        continue
                contained, embedding_count = containment_cache[key]
                if not contained:
                    rejected_field.append((
                        polynomial,
                        "Q[x]/(f) does not embed in H_%s" % D,
                    ))
                    continue
            else:
                embedding_count = None

            orbit_required_bits = holdout_required_bits
            orbit_residuals = [
                _cm_relative_residual_bits(polynomial, value)
                for value in full_numeric_orbit
            ]
            orbit_matches = sum(
                1 for bits in orbit_residuals
                if bits >= orbit_required_bits
            )
            if require_full_orbit and orbit_matches != len(
                full_numeric_orbit
            ):
                rejected_orbit.append((
                    polynomial,
                    orbit_matches,
                    len(full_numeric_orbit),
                    orbit_required_bits,
                ))
                continue

            stable_entry["accepted_windows"].append({
                "window": window,
                "request_degrees": common_request_degrees,
                "orbit_residual_bits": min(orbit_residuals),
            })
            # Keep a fixed algdep degree bound throughout all accepted
            # precisions, and require genuinely separated windows.
            chosen_windows = []
            chosen_degree = None
            for request_degree in sorted(common_request_degrees):
                eligible = [
                    entry for entry in stable_entry["accepted_windows"]
                    if request_degree in entry["request_degrees"]
                ]
                separated = []
                for entry in eligible:
                    center = (entry["window"][0] + entry["window"][-1]) / 2
                    if not separated or all(
                        abs(center - (other["window"][0] + other["window"][-1]) / 2)
                        >= algdep_window_gap_bits
                        for other in separated
                    ):
                        separated.append(entry)
                if len(separated) >= algdep_min_stable_windows:
                    chosen_windows = separated
                    chosen_degree = request_degree
                    break
            if not chosen_windows:
                continue
            stable_precisions = tuple(sorted({
                bits for entry in chosen_windows for bits in entry["window"]
            }))
            accepted_in_window.append({
                "polynomial": polynomial,
                "stable_precisions": stable_precisions,
                "stability_windows": [entry["window"] for entry in chosen_windows],
                "request_degree": chosen_degree,
                "selected_residual_bits": selected_residual_bits,
                "orbit_residual_bits": min(entry["orbit_residual_bits"] for entry in chosen_windows),
                "orbit_matches": orbit_matches,
                "embedding_count": embedding_count,
            })

        if accepted_in_window:
            best = sorted(
                accepted_in_window,
                key=lambda entry: (
                    entry["polynomial"].degree(),
                    _cm_coefficient_height(entry["polynomial"]),
                    -entry["selected_residual_bits"],
                ),
            )[0]
            return {
                "status": "recognized",
                "polynomial": best["polynomial"],
                "known_bits": known_bits,
                "algdep_use_bits": max(best["stable_precisions"]),
                "algdep_stable_precisions": best["stable_precisions"],
                "algdep_stability_windows": best["stability_windows"],
                "algdep_request_degree": best["request_degree"],
                "algdep_holdout_bits": algdep_holdout_bits,
                "algdep_support": len(best["stable_precisions"]),
                "selected_x_residual_bits":
                    best["selected_residual_bits"],
                "orbit_min_residual_bits":
                    best["orbit_residual_bits"],
                "orbit_matches": best["orbit_matches"],
                "orbit_size_with_conjugates":
                    len(full_numeric_orbit),
                "observed_distinct_x": len(
                    _cm_cluster(full_numeric_orbit, known_bits)
                ),
                "class_field_verified":
                    bool(require_class_field_containment),
                "class_field_degree": (
                    ""
                    if class_field is None
                    else Integer(class_field.absolute_degree())
                ),
                "class_field_embedding_count":
                    best["embedding_count"],
                "numerical_warning":
                    selected_item.get("warning", ""),
            }

    attempted_precisions = sorted(candidates_by_precision)
    stable_text = [
        "%s @ %s" % (
            entry["polynomial"],
            entry["windows"],
        )
        for entry in stable_candidates.values()
    ]
    if not stable_candidates:
        status = (
            "unstable_algdep"
            if any(candidates_by_precision.values())
            else "algdep_failed"
        )
        message = (
            "No polynomial was reproduced on any consecutive precision "
            "window. Different algdep answers are rejected."
        )
    elif class_field_error is not None and require_class_field_containment:
        status = "class_field_check_failed"
        message = "Could not construct H_%s: %s" % (D, class_field_error)
    elif rejected_field and not rejected_orbit and not rejected_holdout:
        status = "class_field_mismatch"
        message = (
            "Stable algdep candidates do not define subfields of H_%s. "
            "Rejected: %s" % (
                D,
                "; ".join(
                    "%s (%s)" % (polynomial, reason)
                    for polynomial, reason in rejected_field
                ),
            )
        )
    elif rejected_orbit:
        status = "orbit_mismatch"
        message = (
            "Stable candidate(s) failed the full numerical Galois-orbit "
            "test. Rejected: %s"
            % "; ".join(
                "%s (%s/%s roots at >=%s bits)"
                % values for values in rejected_orbit
            )
        )
    elif rejected_holdout:
        status = "holdout_failed"
        message = (
            "Consecutive-bit candidates failed the %s-bit holdout test "
            "on the original numerical values. No candidate also passed "
            "%s separated stability windows."
            % (algdep_holdout_bits, algdep_min_stable_windows)
        )
    else:
        status = "unstable_algdep"
        message = (
            "No candidate passed %s consecutive-precision windows with "
            "centers separated by at least %s bits at one fixed algdep "
            "degree bound."
            % (algdep_min_stable_windows, algdep_window_gap_bits)
        )

    if errors:
        message += " Last algdep error: %s" % errors[-1]
    return {
        "status": status,
        "polynomial": None,
        "known_bits": known_bits,
        "attempted_use_bits": attempted_precisions,
        "stable_candidates": stable_text,
        "class_field_verified": False,
        "message": message,
    }


def _cm_two_torsion_annotation(E, R, f):
    if f is None:
        return ""
    x = R.gen()
    division_polynomial = (
        4 * x ** 3
        + ZZ(E.b2()) * x ** 2
        + 2 * ZZ(E.b4()) * x
        + ZZ(E.b6())
    )
    return (
        "2-torsion"
        if division_polynomial.change_ring(QQ) % f.change_ring(QQ) == 0
        else ""
    )


def _cm_absolute_degree(field):
    if field is QQ:
        return Integer(1)
    try:
        return Integer(field.absolute_degree())
    except AttributeError:
        return Integer(field.degree())


def _cm_field_polynomial(field):
    if field is QQ:
        return None
    try:
        return field.defining_polynomial()
    except AttributeError:
        return field.polynomial()


def _cm_subfield_description(H, elements, name):
    """
    Return the exact subfield generated by ``elements`` inside H.

    Sage returns the descended elements and the inclusion map into H, so the
    equalities with the original ambient elements can be checked exactly.
    """
    field, descended, inclusion = H.subfield_from_elements(
        list(elements),
        name=name,
        polred=True,
    )
    if len(descended) != len(elements):
        raise ArithmeticError(
            "Subfield descent returned the wrong number of elements."
        )
    if any(
        inclusion(descended[index]) != H(elements[index])
        for index in range(len(elements))
    ):
        raise ArithmeticError(
            "A descended coordinate does not map back to its element in H_D."
        )
    return {
        "field": field,
        "degree": _cm_absolute_degree(field),
        "defining_polynomial": _cm_field_polynomial(field),
        "elements": list(descended),
        "inclusion": inclusion,
    }


def _cm_safe_field_name(prefix, D, operator, index):
    text = "%s_%s_%s_%s" % (
        prefix,
        "m%s" % abs(Integer(D)) if Integer(D) < 0 else Integer(D),
        _cm_sanitize_name(operator),
        Integer(index),
    )
    return re.sub(r"[^A-Za-z0-9_]", "_", str(text))


def _cm_recover_row_points(row, H):
    """
    Recover all exact points associated with one recognized x-polynomial.

    The output contains both signs of y and all Galois-conjugate x-roots.
    """
    E = EllipticCurve(str(row["E"]))
    EH = E.change_ring(H)
    polynomial = row.get("polynomial")
    if polynomial is None:
        raise ValueError("A recognized row has no x-polynomial.")

    RX = PolynomialRing(H, "X")
    x_roots = RX(polynomial).roots(multiplicities=False)
    expected_x_roots = Integer(polynomial.degree())
    if len(x_roots) != expected_x_roots:
        raise ValueError(
            "x-field error: f has %s roots in H_%s, expected %s."
            % (len(x_roots), row["D"], expected_x_roots)
        )

    RY = PolynomialRing(H, "Y")
    Y = RY.gen()
    recovered = []
    point_index = Integer(0)

    for xH in x_roots:
        rhs = (
            xH ** 3
            + H(E.a2()) * xH ** 2
            + H(E.a4()) * xH
            + H(E.a6())
        )
        linear_y = H(E.a1()) * xH + H(E.a3())
        y_polynomial = Y ** 2 + linear_y * Y - rhs
        y_roots = y_polynomial.roots(multiplicities=False)
        if not y_roots:
            raise ValueError(
                "y-field error: the Weierstrass equation has no y-root "
                "in H_%s above x=%s." % (row["D"], xH)
            )

        for yH in y_roots:
            point_index += 1
            point_H = EH([xH, yH, H(1)])

            name_x = _cm_safe_field_name(
                "x", row["D"], row["operator"], point_index
            )
            name_y = _cm_safe_field_name(
                "y", row["D"], row["operator"], point_index
            )
            name_point = _cm_safe_field_name(
                "p", row["D"], row["operator"], point_index
            )

            x_field = _cm_subfield_description(H, [xH], name_x)
            y_field = _cm_subfield_description(H, [yH], name_y)
            point_field = _cm_subfield_description(
                H, [xH, yH], name_point
            )
            L = point_field["field"]
            xL, yL = point_field["elements"]
            EL = E.change_ring(L)
            point_L = EL([xL, yL, L(1)])
            inclusion = point_field["inclusion"]

            mapped_point = EH([
                inclusion(point_L[0]),
                inclusion(point_L[1]),
                H(1),
            ])
            if mapped_point != point_H:
                raise ArithmeticError(
                    "Point-field error: the descended point does not map "
                    "back to the point over H_D."
                )
            if _cm_absolute_degree(H) % point_field["degree"] != 0:
                raise ArithmeticError(
                    "Point-field error: Q(x,y) has degree %s, which does "
                    "not divide [H_D:Q]=%s."
                    % (
                        point_field["degree"],
                        _cm_absolute_degree(H),
                    )
                )

            recovered.append({
                "point_index": point_index,
                "ambient_point": point_H,
                "point": point_L,
                "curve": EL,
                "x": xL,
                "y": yL,
                "x_field": x_field,
                "y_field": y_field,
                "point_field": point_field,
            })

    return {
        "E": str(row["E"]),
        "D": Integer(row["D"]),
        "h": Integer(row["h"]),
        "operator": str(row["operator"]),
        "status": "algebraic_points_recovered",
        "ambient_class_field": H,
        "ambient_class_field_degree": _cm_absolute_degree(H),
        "x_polynomial": polynomial,
        "points": recovered,
        "message": "",
    }
