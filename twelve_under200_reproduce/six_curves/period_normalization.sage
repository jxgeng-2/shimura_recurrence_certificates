"""Match all side periods, with an exact lattice-index calculation."""
from sage.all import ComplexField, RealField, EllipticCurve, ZZ, gcd


def six_lattice_coordinates(z, e1, e2):
    area = e1.real()*e2.imag() - e2.real()*e1.imag()
    return ((z.real()*e2.imag()-e2.real()*z.imag())/area,
            (e1.real()*z.imag()-z.real()*e1.imag())/area)


def six_match_periods(periods, E_label="33a1", prec_bits=1000,
                      search_bound=12, min_match_bits=60):
    """Find a scalar taking every side period into the target period lattice.

    The index is the gcd of ALL two-by-two minors of the integral lattice
    coordinates, not the determinant from an arbitrarily selected pair.
    Smallest index is selected before comparing numerical residuals.
    """
    CC, RR = ComplexField(prec_bits), RealField(prec_bits)
    periods = [CC(p) for p in periods]
    E = EllipticCurve(E_label)
    e1, e2 = map(CC, E.period_lattice().basis(prec=prec_bits))
    tol = RR(2)**(-min_match_bits)
    size = max([abs(p) for p in periods] + [RR(1)])
    nonzero = []
    for p in periods:
        if abs(p) < size*tol:
            continue
        if any(min(abs(p-q), abs(p+q)) < size*tol for q in nonzero):
            continue
        nonzero.append(p)
    nonzero.sort(key=abs)
    candidates, seen_mu = [], []
    for i, p1 in enumerate(nonzero):
        for p2 in nonzero[i+1:]:
            if abs((p2/p1).imag()) < RR("1e-10"):
                continue
            ratio = complex(p2/p1)
            c1, c2 = complex(e1), complex(e2)
            area = c1.real*c2.imag - c2.real*c1.imag
            for m in range(-int(search_bound), int(search_bound)+1):
                for n in range(-int(search_bound), int(search_bound)+1):
                    if m == 0 and n == 0:
                        continue
                    approx = ratio*(m*c1+n*c2)
                    a = (approx.real*c2.imag-c2.real*approx.imag)/area
                    b = (c1.real*approx.imag-approx.real*c1.imag)/area
                    if max(abs(a-round(a)), abs(b-round(b))) > 1e-7:
                        continue
                    mu = (m*e1+n*e2)/p1
                    if mu.real() < 0 or (mu.real() == 0 and mu.imag() < 0):
                        mu = -mu
                    if any(abs(mu-v) < tol*max(1, abs(mu)) for v in seen_mu):
                        continue
                    seen_mu.append(mu)
                    coords = [six_lattice_coordinates(mu*p, e1, e2)
                              for p in periods]
                    integers = [(ZZ(a.round()), ZZ(b.round())) for a,b in coords]
                    error = max([max(abs(a-u), abs(b-v))
                                 for (a,b),(u,v) in zip(coords,integers)] + [RR(0)])
                    if error > tol:
                        continue
                    index = gcd([abs(u*s-v*r) for u,v in integers for r,s in integers])
                    if index == 0:
                        continue
                    candidates.append(dict(mu=mu, index=index,
                        lattice_coordinates=integers, lattice_error=error,
                        lattice_error_bits=(-error.log(2)).floor() if error else prec_bits,
                        E_label=E_label, periods=periods, target_periods=(e1,e2)))
            primitive = [v for v in candidates if v["index"] == 1]
            if primitive:
                return min(primitive, key=lambda v: v["lattice_error"])
    if not candidates:
        raise ValueError("No all-period lattice match at %s bits for %s" %
                         (min_match_bits, E_label))
    return min(candidates, key=lambda v: (v["index"], v["lattice_error"]))
