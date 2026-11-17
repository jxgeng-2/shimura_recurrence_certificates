"""Reproducible N=21, D=-43 CM series and exact side-generator geometry.

The recurrence and scale are the inputs in cm_expansions_voight_willis.md.
The attached notebook independently checks the initial values/modularity to
about 28 decimal places; high precision values below are generated afresh.
"""

from sage.all import QQ, ZZ, RealField, ComplexField, matrix, kronecker
import json


def x21_quaternion_matrix(coords, RF):
    x, y, z, t = [QQ(c) for c in coords]
    s = RF(21).sqrt()
    return matrix(RF, [[RF(x)+RF(y)*s, RF(z)+RF(t)*s],
                       [-RF(z)+RF(t)*s, RF(x)-RF(y)*s]])


def x21_cm43_tau(CC):
    RF = RealField(CC.precision())
    return CC((-RF(21)-5*RF(21).sqrt())/16,
              (5*RF(43).sqrt()+RF(903).sqrt())/16)


def x21_cm43_theta(CC):
    """Theta = pi Omega_-43^2 (-63/4 - 3 sqrt(-903)/4)."""
    RF = RealField(CC.precision())
    log_product = sum(ZZ(kronecker(-43, j))*(RF(j)/43).log_gamma()
                      for j in range(1, 43))
    omega_squared = log_product.exp()/(2*RF.pi()*43)
    alpha = CC(-RF(63)/4, -3*RF(903).sqrt()/4)
    return CC(RF.pi()*omega_squared)*alpha


def x21_cm43_series(Nper=1000, prec_bits=1000, guard_bits=None):
    """Generate b[0:Nper] by a quadratic-time coefficient recurrence.

    Coefficients in the real q=Theta*w coordinate are computed with guard
    precision: the auxiliary elliptic coordinate is meromorphic and can
    magnify intermediate roundoff. Compare two guarded runs if increasing
    Nper substantially. The returned coefficients have prec_bits precision.
    """
    Nper, prec_bits = int(Nper), int(prec_bits)
    if Nper < 2:
        raise ValueError('Nper must be at least 2.')
    if guard_bits is None:
        guard_bits = 3*Nper + 128
    work_bits = prec_bits + int(guard_bits)
    RF = RealField(work_bits)
    CC = ComplexField(work_bits)
    out_CC = ComplexField(prec_bits)
    lam = -RF(43)/3969
    u0, u1 = RF(441)/25, -RF(9632)/375

    # A'=B, U'=V. C=1/A, T=1/(U-49). All arrays are ordinary
    # (not exponential) Taylor coefficients in q.
    A, B, U, V = [RF(1)], [RF(1)], [u0], [u1]
    C, T = [RF(1)], [1/(u0-49)]
    A2, A3, B2, BV, U2, T2, numerator, rpart = ([] for _ in range(8))

    def conv(a, b, n):
        return sum((a[k]*b[n-k] for k in range(n+1)), RF(0))

    for n in range(Nper-1):
        A2.append(conv(A, A, n))
        A3.append(conv(A2, A, n))
        B2.append(conv(B, B, n))
        BV.append(conv(B, V, n))
        U2.append(conv(U, U, n))
        T2.append(conv(T, T, n))
        numerator.append(17*U2[n]+392*U[n]-(2401 if n == 0 else 0))
        rpart.append(conv(numerator, T2, n))
        rhsA = RF(3)/2*conv(B2, C, n)+16*lam*conv(rpart, A3, n)
        # H'(U) = -84 U^2 +752 U -1372.
        hprime_a2 = -84*conv(U2, A2, n)+752*conv(U, A2, n)-1372*A2[n]
        rhsU = conv(BV, C, n)+lam/2*hprime_a2
        A.append(B[n]/(n+1))
        U.append(V[n]/(n+1))
        B.append(rhsA/(n+1))
        V.append(rhsU/(n+1))
        C.append(-sum((A[k]*C[n+1-k] for k in range(1,n+2)), RF(0)))
        T.append(-sum((U[k]*T[n+1-k] for k in range(1,n+2)), RF(0))/(u0-49))

    theta = x21_cm43_theta(CC)
    theta_power = CC(1)
    b = []
    for an in A:
        b.append(out_CC(an*theta_power))
        theta_power *= theta
    return {'b': b, 'theta': out_CC(theta), 'working_bits': work_bits,
            'prec_bits': prec_bits, 'Nper': Nper,
            'guard_bits': int(guard_bits)}


def build_x21_cm43_setup(side_data_path, Nper=1000, prec_bits=1000,
                        guard_bits=None):
    """Build series, tau0, side matrices directly from exact JSON entries."""
    with open(side_data_path, encoding='utf-8-sig') as stream:
        data = json.load(stream)
    if (QQ(data['quaternion_a']), QQ(data['quaternion_b'])) != (21,-1):
        raise ValueError('Expected quaternion presentation (21,-1).')
    if QQ(data['cm_discriminant']) != -43:
        raise ValueError('Expected CM center discriminant -43.')
    RF, CC = RealField(prec_bits), ComplexField(prec_bits)
    out = x21_cm43_series(Nper, prec_bits, guard_bits)
    out['tau0'] = x21_cm43_tau(CC)
    out['S_matrices'] = [x21_quaternion_matrix(c, RF)
                         for c in data['side_generators']]
    out['side_data_path'] = str(side_data_path)
    out['side_data'] = data
    out['E_label'] = '21a1'
    return out


def audit_x21_series_modularity(setup, Nper=None):
    """Test every side at the hyperbolic midpoint of tau and g^-1(tau).

    Both evaluation points have the same disk radius, minimizing the worst
    series radius along that pair. Residuals test the weight-two law.
    """
    b, tau = setup['b'], setup['tau0']
    CC = tau.parent()
    RF = RealField(CC.precision())
    if Nper is None:
        Nper = len(b)
    Nper = min(int(Nper), len(b))
    tbar = tau.conjugate()

    def mob(g, z):
        return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])

    def disk(z):
        return (z-tau)/(z-tbar)

    def evaluate(w):
        value = CC(0)
        for an in reversed(b[:Nper]):
            value = value*w+an
        return (1-w)**2*value

    rows = []
    for index, g in enumerate(setup['S_matrices']):
        endpoint = disk(mob(g.inverse(), tau))
        w = endpoint/(1+(1-abs(endpoint)**2).sqrt())
        z = (tau-tbar*w)/(1-w)
        gz = mob(g,z)
        gw = disk(gz)
        left = evaluate(gw)
        right = (g[1,0]*z+g[1,1])**2/g.det()*evaluate(w)
        scale = max(RF(1), abs(left), abs(right))
        error = abs(left-right)/scale
        bits = CC.precision() if error == 0 else -error.log2()
        rows.append({'side_index': index, 'radius': max(abs(w),abs(gw)),
                     'relative_error':error, 'residual_bits':bits})
    return rows
