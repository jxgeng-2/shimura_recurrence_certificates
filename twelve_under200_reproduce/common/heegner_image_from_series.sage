"""CM divisor images from a disk expansion, with reusable Dirichlet reduction.

The side-data JSON is read without modifying it. Exact quaternion coordinates
are preferred to exported decimal matrices. Endpoint reduction computes an
integral modulo source periods; mu must send these periods into the target
elliptic-curve lattice. Use method='path' for actual path integrals/periods.

    reducer = DirichletReducer('x21_cm43_side_data.json', prec_bits=1000)
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
