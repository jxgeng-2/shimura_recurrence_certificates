"""N=15, D=-7 recurrence, exact geometry, and independent modularity tests.

The recurrence is an input from cm_expansions_voight_willis.md, not a
standalone proof of modularity.  Always inspect the side-pairing audit.
"""
from sage.all import QQ, ZZ, RealField, ComplexField, matrix, kronecker
import json


def x15_quaternion_matrix(coords, RF):
    x, y, z, t = [QQ(c) for c in coords]
    s, r = RF(15).sqrt(), RF(30).sqrt()
    return matrix(RF, [[RF(x)+RF(y)*s, (RF(z)+RF(t)*s)*r],
                       [(-RF(z)+RF(t)*s)*r, RF(x)-RF(y)*s]])


def x15_cm7_tau(CC):
    RF = RealField(CC.precision())
    return CC(0, (RF(210).sqrt()-RF(14).sqrt())/14)


def x15_cm7_theta(CC):
    """The exact lift has Theta=-8*sqrt(15)/15*pi*Omega_-7^2."""
    RF=RealField(CC.precision())
    lp=sum(ZZ(kronecker(-7,j))*(RF(j)/7).log_gamma() for j in range(1,7))
    return CC(-RF(8)/15*RF(15).sqrt()*lp.exp()/14)


def x15_cm7_series(Nper=1000, prec_bits=1000, guard_bits=None):
    """Quadratic-time guarded recurrence for ordinary Taylor coefficients."""
    Nper, prec_bits = int(Nper), int(prec_bits)
    if Nper < 2:
        raise ValueError('Nper must be at least 2.')
    if guard_bits is None:
        guard_bits = 3*Nper+128
    work_bits = prec_bits+int(guard_bits)
    RF, CC = RealField(work_bits), ComplexField(work_bits)
    out_CC = ComplexField(prec_bits)
    lam = -RF(7)/64
    u0,u1=RF(9),-RF(63)
    A,B,U,V=[RF(1)],[RF(1)],[u0],[u1]
    C=[RF(1)]
    A2,A3,B2,BV,U2=([] for _ in range(5))
    def conv(a,b,n):
        return sum((a[k]*b[n-k] for k in range(n+1)),RF(0))
    for n in range(Nper-1):
        A2.append(conv(A,A,n))
        A3.append(conv(A2,A,n))
        B2.append(conv(B,B,n))
        BV.append(conv(B,V,n))
        U2.append(conv(U,U,n))
        rhsA=RF(3)/2*conv(B2,C,n)+lam*(4*A3[n]+RF(4)/3*conv(U,A3,n))
        hp=-36*conv(U2,A2,n)-656*conv(U,A2,n)-108*A2[n]
        rhsU=conv(BV,C,n)+lam/2*hp
        A.append(B[n]/(n+1))
        U.append(V[n]/(n+1))
        B.append(rhsA/(n+1))
        V.append(rhsU/(n+1))
        C.append(-sum((A[k]*C[n+1-k] for k in range(1,n+2)),RF(0)))
    theta = x15_cm7_theta(CC)
    power = CC(1)
    b = []
    for an in A:
        b.append(out_CC(an*power))
        power *= theta
    return dict(b=b,theta=out_CC(theta),working_bits=work_bits,
                prec_bits=prec_bits,Nper=Nper,guard_bits=int(guard_bits))


def build_x15_cm7_setup(side_data_path, Nper=1000, prec_bits=1000,
                         guard_bits=None):
    with open(side_data_path,encoding='utf-8-sig') as stream:
        data=json.load(stream)
    if (QQ(data['quaternion_a']),QQ(data['quaternion_b'])) != (15,-30):
        raise ValueError('Expected quaternion presentation (15,-30).')
    if QQ(data['cm_discriminant']) != -7:
        raise ValueError('Expected CM center discriminant -7.')
    if list(map(QQ,data['cm_element'])) != [QQ(1)/2,0,QQ(1)/4,-QQ(1)/60]:
        raise ValueError('The supplied recurrence is for the specified -7 lift.')
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out=x15_cm7_series(Nper,prec_bits,guard_bits)
    out.update(tau0=x15_cm7_tau(CC),
               S_matrices=[x15_quaternion_matrix(c,RF) for c in data['side_generators']],
               side_data_path=str(side_data_path),side_data=data,E_label='15a1')
    return out


def audit_x15_series_modularity(setup,Nper=None):
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


