"""N=46, D=-163 recurrence, exact geometry, and independent modularity tests.

The recurrence is an input from cm_expansions_voight_willis.md, not a
standalone proof of modularity.  Always inspect the side-pairing audit.
"""
from sage.all import QQ, ZZ, RealField, ComplexField, matrix, kronecker
import json


def x46_quaternion_matrix(coords, RF):
    x, y, z, t = [QQ(c) for c in coords]
    s, r = RF(46).sqrt(), RF(3).sqrt()
    return matrix(RF, [[RF(x)+RF(y)*s, (RF(z)+RF(t)*s)*r],
                       [(-RF(z)+RF(t)*s)*r, RF(x)-RF(y)*s]])


def x46_cm163_tau(CC):
    RF = RealField(CC.precision())
    return CC(0, (35*RF(489).sqrt()+4*RF(22494).sqrt())/489)


def x46_cm163_theta(CC):
    """Theta = pi Omega_-163^2 * (7855/69)*sqrt(46)."""
    RF = RealField(CC.precision())
    log_product = sum(ZZ(kronecker(-163,j))*(RF(j)/163).log_gamma()
                      for j in range(1,163))
    omega_squared = log_product.exp()/(2*RF.pi()*163)
    return CC(RF.pi()*omega_squared*(RF(7855)/69)*RF(46).sqrt())


def x46_cm163_series(Nper=1000, prec_bits=1000, guard_bits=None):
    """Quadratic-time guarded recurrence for ordinary Taylor coefficients."""
    Nper, prec_bits = int(Nper), int(prec_bits)
    if Nper < 2:
        raise ValueError('Nper must be at least 2.')
    if guard_bits is None:
        guard_bits = 3*Nper+128
    work_bits = prec_bits+int(guard_bits)
    RF, CC = RealField(work_bits), ComplexField(work_bits)
    out_CC = ComplexField(prec_bits)
    lam = -RF(844992)/61701025
    u0, u1 = RF(3136)/81, -RF(15188992)/636255
    A, B, U, V = [RF(1)], [RF(1)], [u0], [u1]
    C, T = [RF(1)], [1/(u0-64)]
    A2,A3,B2,BV,U2,U3,T2,numerator,rpart = ([] for _ in range(9))

    def conv(a,b,n):
        return sum((a[k]*b[n-k] for k in range(n+1)),RF(0))

    for n in range(Nper-1):
        A2.append(conv(A,A,n))
        A3.append(conv(A2,A,n))
        B2.append(conv(B,B,n))
        BV.append(conv(B,V,n))
        U2.append(conv(U,U,n))
        U3.append(conv(U2,U,n))
        T2.append(conv(T,T,n))
        numerator.append(3*U3[n]+52*U2[n]+13824*U[n]-(311296 if n==0 else 0))
        rpart.append(conv(numerator,T2,n))
        rhsA = RF(3)/2*conv(B2,C,n)+lam/8*conv(rpart,A3,n)
        # H(u)=-4u^3+180u^2-2048u, so H'=-12u^2+360u-2048.
        hprime_a2 = -12*conv(U2,A2,n)+360*conv(U,A2,n)-2048*A2[n]
        rhsU = conv(BV,C,n)+lam/2*hprime_a2
        A.append(B[n]/(n+1))
        U.append(V[n]/(n+1))
        B.append(rhsA/(n+1))
        V.append(rhsU/(n+1))
        C.append(-sum((A[k]*C[n+1-k] for k in range(1,n+2)),RF(0)))
        T.append(-sum((U[k]*T[n+1-k] for k in range(1,n+2)),RF(0))/(u0-64))
    theta = x46_cm163_theta(CC)
    power = CC(1)
    b = []
    for an in A:
        b.append(out_CC(an*power))
        power *= theta
    return dict(b=b,theta=out_CC(theta),working_bits=work_bits,
                prec_bits=prec_bits,Nper=Nper,guard_bits=int(guard_bits))


def build_x46_cm163_setup(side_data_path, Nper=1000, prec_bits=1000,
                         guard_bits=None):
    with open(side_data_path,encoding='utf-8-sig') as stream:
        data=json.load(stream)
    if (QQ(data['quaternion_a']),QQ(data['quaternion_b'])) != (46,-3):
        raise ValueError('Expected quaternion presentation (46,-3).')
    if QQ(data['cm_discriminant']) != -163:
        raise ValueError('Expected CM center discriminant -163.')
    if list(map(QQ,data['cm_element'])) != [QQ(1)/2,0,-QQ(35)/6,-QQ(2)/3]:
        raise ValueError('The supplied recurrence is for the specified -163 lift.')
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out=x46_cm163_series(Nper,prec_bits,guard_bits)
    out.update(tau0=x46_cm163_tau(CC),
               S_matrices=[x46_quaternion_matrix(c,RF) for c in data['side_generators']],
               side_data_path=str(side_data_path),side_data=data,E_label='46a1')
    return out


def audit_x46_series_modularity(setup,Nper=None):
    """Weight-two tests at balanced disk points for every exact side matrix."""
    b,tau=setup['b'],setup['tau0']
    CC=tau.parent()
    RF=RealField(CC.precision())
    Nper=min(len(b),len(b) if Nper is None else int(Nper))
    tbar=tau.conjugate()
    def mob(g,z):
        return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])
    def disk(z):
        return (z-tau)/(z-tbar)
    def evaluate(w):
        value=CC(0)
        for an in reversed(b[:Nper]):
            value=value*w+an
        return (1-w)**2*value
    rows=[]
    for index,g in enumerate(setup['S_matrices']):
        endpoint=disk(mob(g.inverse(),tau))
        w=endpoint/(1+(1-abs(endpoint)**2).sqrt())
        z=(tau-tbar*w)/(1-w)
        gz=mob(g,z)
        gw=disk(gz)
        left=evaluate(gw)
        right=(g[1,0]*z+g[1,1])**2/g.det()*evaluate(w)
        error=abs(left-right)/max(RF(1),abs(left),abs(right))
        bits=CC.precision() if error==0 else -error.log2()
        rows.append(dict(side_index=index,radius=max(abs(w),abs(gw)),
                         relative_error=error,residual_bits=bits))
    return rows
