"""Guarded disk series for all twelve forms; exact inputs stay inline."""
# six_curves\x14_cm11_setup.sage
"""Guarded X(14,1) disk expansion at the attached CM lift of D=-11."""
from sage.all import *
import json


def x14_quaternion_matrix(coords,RF):
    x,y,z,t=map(QQ,coords)
    s,r=RF(7).sqrt(),RF(2).sqrt()
    return matrix(RF,[[x+y*s,(z+t*s)*r],[(-z+t*s)*r,x-y*s]])


def x14_cm11_tau(CC):
    RF=RealField(CC.precision())
    return CC(-RF(14).sqrt()/6,RF(22).sqrt()/6)


def x14_cm11_theta(CC):
    RF=RealField(CC.precision())
    log_product=sum(ZZ(kronecker(-11,j))*(RF(j)/11).log_gamma() for j in range(1,11))
    omega2=log_product.exp()/(2*RF.pi()*11)
    return CC(RF.pi()*omega2)*CC(1,RF(77).sqrt()/7)


def x14_cm11_series(Nper=1400,prec_bits=1000,guard_bits=None):
    Nper,prec_bits=int(Nper),int(prec_bits)
    guard_bits=3*Nper+128 if guard_bits is None else int(guard_bits)
    work_bits=prec_bits+guard_bits
    RF,CC=RealField(work_bits),ComplexField(work_bits)
    lam=-RF(11)/9
    A,B,U,V,C=[RF(1)],[RF(1)],[RF(16)],[-RF(352)/3],[RF(1)]
    A2,A3,B2,BV,U2=([] for _ in range(5))
    def conv(a,b,n):
        return sum((a[k]*b[n-k] for k in range(n+1)),RF(0))
    for n in range(Nper-1):
        A2.append(conv(A,A,n)); A3.append(conv(A2,A,n))
        B2.append(conv(B,B,n)); BV.append(conv(B,V,n)); U2.append(conv(U,U,n))
        rhsA=RF(3)/2*conv(B2,C,n)-lam*RF(3)/8*(4*A3[n]-conv(U,A3,n))
        rhsU=conv(BV,C,n)+lam/2*(-12*conv(U2,A2,n)+104*conv(U,A2,n)-512*A2[n])
        A.append(B[n]/(n+1)); U.append(V[n]/(n+1))
        B.append(rhsA/(n+1)); V.append(rhsU/(n+1))
        C.append(-sum((A[k]*C[n+1-k] for k in range(1,n+2)),RF(0)))
    theta=x14_cm11_theta(CC)
    outCC=ComplexField(prec_bits)
    power=CC(1); coeffs=[]
    for v in A:
        coeffs.append(outCC(v*power));power*=theta
    return dict(b=coeffs,theta=outCC(theta),working_bits=work_bits,
                prec_bits=prec_bits,Nper=Nper,guard_bits=guard_bits)


def build_x14_cm11_setup(side_data_path=None,Nper=1400,prec_bits=1000,guard_bits=None):
    if side_data_path is None: side_data_path = side_data(14)
    data = _heegner_read_side_data(side_data_path)
    assert (QQ(data['quaternion_a']),QQ(data['quaternion_b']))==(7,-2)
    assert list(map(QQ,data['cm_element']))==[QQ(1)/2,-QQ(1)/2,-QQ(3)/2,0]
    out=x14_cm11_series(Nper,prec_bits,guard_bits)
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out.update(tau0=x14_cm11_tau(CC),S_matrices=[x14_quaternion_matrix(c,RF) for c in data['side_generators']],
               side_data=data,side_data_path=str(side_data_path),E_label='14a1')
    return out


def audit_x14_series_modularity(setup,Nper=None):
    coeffs,tau=setup['b'],setup['tau0']; CC=tau.parent(); RF=RealField(CC.precision())
    nper=len(coeffs) if Nper is None else min(int(Nper),len(coeffs));bar=tau.conjugate()
    def mob(g,z): return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])
    def disk(z): return (z-tau)/(z-bar)
    def val(w):
        s=CC(0)
        for b in reversed(coeffs[:nper]):s=s*w+b
        return (1-w)**2*s
    rows=[]
    for i,g in enumerate(setup['S_matrices']):
        endpoint=disk(mob(g.inverse(),tau))
        w=endpoint/(1+(1-abs(endpoint)**2).sqrt());z=(tau-bar*w)/(1-w)
        gw=disk(mob(g,z));left=val(gw);right=(g[1,0]*z+g[1,1])**2/g.det()*val(w)
        err=abs(left-right)/max(RF(1),abs(left),abs(right))
        rows.append(dict(side_index=i,radius=max(abs(w),abs(gw)),relative_error=err,
                         residual_bits=CC.precision() if err==0 else -err.log2()))
    return rows


# six_curves\x15_cm7_setup.sage
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
    data = _heegner_read_side_data(side_data_path)
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




# six_curves\x21_cm43_setup.sage
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
    data = _heegner_read_side_data(side_data_path)
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


# genus_two\x26_cm11_setup.sage
"""Guarded ordinary-Taylor recurrences from the Tables 2 and 7.

The six-dimensional first-order system retains the elliptic quotient
coordinate separately for each eigenform. No scalar closure is assumed.
"""
from sage.all import *
from pathlib import Path
import time

X26_GAMMA_COORDS = [[5,-1,QQ(6)/5,QQ(1)/5],[5,-1,-QQ(6)/5,-QQ(1)/5],
    [QQ(5)/2,-QQ(1)/2,-QQ(1)/2,0],[QQ(3)/2,0,QQ(1)/10,QQ(1)/10],
    [5,1,-QQ(6)/5,-QQ(1)/5],[5,1,QQ(6)/5,QQ(1)/5],
    [QQ(5)/2,QQ(1)/2,QQ(1)/2,0],[QQ(3)/2,0,-QQ(1)/10,-QQ(1)/10]]

def x26_quaternion_matrix(coords, RF):
    return _heegner_quaternion_matrix(coords,26,-5,RF)

def x26_cm11_tau(CC):
    return _heegner_fixed_point(x26_quaternion_matrix([QQ(1)/2,1,-QQ(7)/2,-QQ(1)/2],RealField(CC.precision())),CC)

def x26_cm11_theta(form,CC):
    RF=RealField(CC.precision())
    lp=sum(ZZ(kronecker(-11,j))*(RF(j)/11).log_gamma() for j in range(1,11))
    omega2=lp.exp()/(2*RF.pi()*11)
    alpha=CC(988-33*RF(26).sqrt(),-156*RF(11).sqrt()-19*RF(286).sqrt())/2990
    return (-21 if form=='26a' else 1)*RF.pi()*omega2*alpha

def x26_cm11_series(form='26a',Nper=1000,prec_bits=1000,guard_bits=None):
    assert form in ('26a','26b')
    Nper,prec_bits=int(Nper),int(prec_bits)
    if guard_bits is None:guard_bits=4*Nper+128
    wb=prec_bits+int(guard_bits);RF,CC=RealField(wb),ComplexField(wb)
    if form=='26a':
        L=-RF(11)/1764;u1=RF(44)/21;r0,r1,r2=RF(9)/2,RF(-9),-RF(507)/2
        hp0,hp1,hp2=map(RF,[-96,152,-24]);known=[1,1,-QQ(16)/147,QQ(274)/3087,-QQ(2488)/9261,QQ(352768)/1361367,QQ(5961776)/28588707]
    else:
        L=-RF(11)/4;u1=RF(44);r0,r1,r2=map(RF,[3,14,-3])
        hp0,hp1,hp2=map(RF,[76,-192,-2028]);known=[1,1,40,-734,-11240,1308544,5134480]
    A,p,W,U,V,J=[[RF(1)],[RF(1)],[RF(1)],[RF(1)],[u1],[RF(1)]]
    J2,R,Hprime=([] for _ in range(3))
    def conv(a,b,n):return sum((a[k]*b[n-k] for k in range(n+1)),RF(0))
    for n in range(Nper-1):
        J2.append(conv(J,J,n))
        R.append((r0 if n==0 else RF(0))+r1*J[n]+r2*J2[n])
        Hprime.append((hp0 if n==0 else RF(0))+hp1*U[n]+hp2*conv(U,U,n))
        nextA=conv(p,A,n)/(n+1)
        nextp=(conv(p,p,n)/2-L*conv(R,W,n))/(n+1)
        nextW=2*conv(p,W,n)/(n+1)
        nextU=V[n]/(n+1)
        nextV=(conv(p,V,n)+L/2*conv(Hprime,W,n))/(n+1)
        nextJ=-conv(V,J2,n)/(n+1)
        for seq,v in zip((A,p,W,U,V,J),(nextA,nextp,nextW,nextU,nextV,nextJ)):seq.append(v)
    for n,c in enumerate(known[:Nper]):
        assert abs(A[n]*factorial(n)-RF(c)) < RF(2)**(-prec_bits)*max(1,abs(RF(c)))
    theta=x26_cm11_theta(form,CC);power=CC(1);out=ComplexField(prec_bits);b=[]
    for an in A:b.append(out(an*power));power*=theta
    return dict(b=b,theta=out(theta),form=form,Nper=Nper,prec_bits=prec_bits,guard_bits=int(guard_bits),working_bits=wb)

def build_x26_cm11_setup(form='26a',Nper=1000,prec_bits=1000,guard_bits=None):
    side=side_data(26)
    assert (QQ(side['quaternion_a']),QQ(side['quaternion_b']))==(26,-5)
    assert list(map(QQ,side['cm_element']))==[QQ(1)/2,1,-QQ(7)/2,-QQ(1)/2]
    out=x26_cm11_series(form,Nper,prec_bits,guard_bits)
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out.update(tau0=x26_cm11_tau(CC),S_matrices=[x26_quaternion_matrix(c,RF) for c in side['side_generators']],
        gamma_matrices=[x26_quaternion_matrix(c,RF) for c in X26_GAMMA_COORDS],
        side_data=side,side_data_path=None,E_label=form+'1')
    return out

def audit_x26_series_modularity(setup,Nper=None,gamma_only=False):
    b,tau=setup['b'],setup['tau0'];CC=tau.parent();RF=RealField(CC.precision());tbar=tau.conjugate()
    Nper=min(len(b),len(b) if Nper is None else int(Nper))
    def mob(g,z):return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])
    def disk(z):return (z-tau)/(z-tbar)
    def evaluate(w):
        value=CC(0)
        for an in reversed(b[:Nper]):value=value*w+an
        return (1-w)**2*value
    rows=[]
    matrices=setup['gamma_matrices'] if gamma_only else setup['S_matrices']
    for i,g in enumerate(matrices):
        end=disk(mob(g.inverse(),tau));w=end/(1+(1-abs(end)**2).sqrt())
        z=(tau-tbar*w)/(1-w);gz=mob(g,z);gw=disk(gz)
        left=evaluate(gw);right=(g[1,0]*z+g[1,1])**2/g.det()*evaluate(w)
        error=abs(left-right)/max(RF(1),abs(left),abs(right))
        rows.append(dict(side_index=i,radius=max(abs(w),abs(gw)),relative_error=error,residual_bits=CC.precision() if error==0 else -error.log2()))
    return rows


# six_curves\x33_cm67_setup.sage
"""Guarded X(33,1) series at the exact D=-67 lift; audit before use."""
from sage.all import QQ,ZZ,RealField,ComplexField,matrix,kronecker
import json

def x33_quaternion_matrix(coords,RF):
    x,y,z,t=map(QQ,coords)
    s=RF(33).sqrt()
    return matrix(RF,[[RF(x)+RF(y)*s,RF(z)+RF(t)*s],
                       [-RF(z)+RF(t)*s,RF(x)-RF(y)*s]])

def x33_cm67_tau(CC):
    RF=RealField(CC.precision())
    return CC(-RF(33).sqrt()/10,RF(67).sqrt()/10)

def x33_cm67_theta(CC):
    RF=RealField(CC.precision())
    lp=sum(ZZ(kronecker(-67,j))*(RF(j)/67).log_gamma() for j in range(1,67))
    omega2=lp.exp()/(2*RF.pi()*67)
    alpha=CC(-RF(23)/5,-RF(23)*RF(2211).sqrt()/165)
    return CC(RF.pi()*omega2)*alpha

def x33_cm67_series(Nper=1000,prec_bits=1000,guard_bits=None):
    Nper,prec_bits=int(Nper),int(prec_bits)
    if guard_bits is None: guard_bits=4*Nper+128
    wb=prec_bits+int(guard_bits)
    RF,CC=RealField(wb),ComplexField(wb)
    lam=-RF(67)/529
    u0,u1=RF(9)/25,-RF(19296)/2875
    A,B,U,V=[RF(1)],[RF(1)],[u0],[u1]
    C,S,T=[RF(1)],[1/u0],[1/(u0-1)]
    A2,A3,B2,BV,U2,T2,rpart=([] for _ in range(7))
    def conv(a,b,n):
        return sum((a[k]*b[n-k] for k in range(n+1)),RF(0))
    for n in range(Nper-1):
        A2.append(conv(A,A,n)); A3.append(conv(A2,A,n))
        B2.append(conv(B,B,n)); BV.append(conv(B,V,n))
        U2.append(conv(U,U,n)); T2.append(conv(T,T,n))
        # R(u)=-4*(1+27/u+24/(u-1)+96/(u-1)^2).
        rpart.append((1 if n==0 else 0)+27*S[n]+24*T[n]+96*T2[n])
        rhsA=RF(3)/2*conv(B2,C,n)+4*lam*conv(rpart,A3,n)
        rhsU=conv(BV,C,n)+lam/2*(-36*conv(U2,A2,n)-80*conv(U,A2,n)-972*A2[n])
        A.append(B[n]/(n+1)); U.append(V[n]/(n+1))
        B.append(rhsA/(n+1)); V.append(rhsU/(n+1))
        C.append(-sum((A[k]*C[n+1-k] for k in range(1,n+2)),RF(0)))
        S.append(-sum((U[k]*S[n+1-k] for k in range(1,n+2)),RF(0))/u0)
        T.append(-sum((U[k]*T[n+1-k] for k in range(1,n+2)),RF(0))/(u0-1))
    theta=x33_cm67_theta(CC)
    power=CC(1); out=ComplexField(prec_bits); b=[]
    for an in A:
        b.append(out(an*power)); power*=theta
    return dict(b=b,theta=out(theta),Nper=Nper,prec_bits=prec_bits,
                guard_bits=int(guard_bits),working_bits=wb)

def build_x33_cm67_setup(side_data_path=None,Nper=1000,prec_bits=1000,guard_bits=None):
    if side_data_path is None: side_data_path = side_data(33)
    data = _heegner_read_side_data(side_data_path)
    assert (QQ(data['quaternion_a']),QQ(data['quaternion_b']))==(33,-1)
    assert list(map(QQ,data['cm_element']))==[QQ(1)/2,-QQ(1)/2,-5,0]
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out=x33_cm67_series(Nper,prec_bits,guard_bits)
    out.update(tau0=x33_cm67_tau(CC),S_matrices=[x33_quaternion_matrix(c,RF) for c in data['side_generators']],
               side_data=data,side_data_path=str(side_data_path),E_label='33a1')
    return out

def audit_x33_series_modularity(setup,Nper=None):
    b,tau=setup['b'],setup['tau0']; CC=tau.parent(); RF=RealField(CC.precision())
    Nper=min(len(b),len(b) if Nper is None else int(Nper)); tbar=tau.conjugate()
    def mob(g,z): return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])
    def disk(z): return (z-tau)/(z-tbar)
    def evaluate(w):
        value=CC(0)
        for an in reversed(b[:Nper]): value=value*w+an
        return (1-w)**2*value
    rows=[]
    for index,g in enumerate(setup['S_matrices']):
        endpoint=disk(mob(g.inverse(),tau)); w=endpoint/(1+(1-abs(endpoint)**2).sqrt())
        z=(tau-tbar*w)/(1-w); gz=mob(g,z); gw=disk(gz)
        left=evaluate(gw); right=(g[1,0]*z+g[1,1])**2/g.det()*evaluate(w)
        error=abs(left-right)/max(RF(1),abs(left),abs(right))
        bits=CC.precision() if error==0 else -error.log2()
        rows.append(dict(side_index=index,radius=max(abs(w),abs(gw)),relative_error=error,residual_bits=bits))
    return rows


# six_curves\x34_cm11_setup.sage
"""N=34, D=-11 recurrence, exact geometry, and independent modularity tests.

The recurrence is an input from cm_expansions_voight_willis.md, not a
standalone proof of modularity.  Always inspect the side-pairing audit.
"""
from sage.all import QQ, ZZ, RealField, ComplexField, matrix, kronecker
import json


def x34_quaternion_matrix(coords, RF):
    x, y, z, t = [QQ(c) for c in coords]
    s, r = RF(34).sqrt(), RF(3).sqrt()
    return matrix(RF, [[RF(x)+RF(y)*s, (RF(z)+RF(t)*s)*r],
                       [(-RF(z)+RF(t)*s)*r, RF(x)-RF(y)*s]])


def x34_cm11_tau(CC):
    RF = RealField(CC.precision())
    return CC((-34*RF(3).sqrt()+7*RF(102).sqrt())/45,
              (7*RF(33).sqrt()-RF(1122).sqrt())/45)


def x34_cm11_theta(CC):
    """Theta with the corrected signed N34 factor in the source note."""
    RF = RealField(CC.precision())
    log_product = sum(ZZ(kronecker(-11,j))*(RF(j)/11).log_gamma()
                      for j in range(1,11))
    omega_squared = log_product.exp()/(2*RF.pi()*11)
    alpha = CC(1-11*RF(34).sqrt()/102, -RF(11).sqrt()/3-RF(374).sqrt()/34)
    return CC(RF.pi()*omega_squared)*alpha


def x34_cm11_series(Nper=1000, prec_bits=1000, guard_bits=None):
    """Quadratic-time guarded recurrence for ordinary Taylor coefficients."""
    Nper, prec_bits = int(Nper), int(prec_bits)
    if Nper < 2:
        raise ValueError('Nper must be at least 2.')
    if guard_bits is None:
        guard_bits = 3*Nper+128
    work_bits = prec_bits+int(guard_bits)
    RF, CC = RealField(work_bits), ComplexField(work_bits)
    out_CC = ComplexField(prec_bits)
    lam = -RF(11)/25
    u0, u1 = -RF(15)/4, -RF(3179)/80
    A, B, U, V = [RF(1)], [RF(1)], [u0], [u1]
    C, T = [RF(1)], [1/(15*u0-16)]
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
        numerator.append(279*U2[n]+3104*U[n]-(39264 if n==0 else 0))
        rpart.append(conv(numerator,T2,n))
        rhsA = RF(3)/2*conv(B2,C,n)-lam/2*conv(rpart,A3,n)
        # H=(u^2+4)(-3u^2+26u-59).
        hprime_a2 = -12*conv(U3,A2,n)+78*conv(U2,A2,n)-142*conv(U,A2,n)+104*A2[n]
        rhsU = conv(BV,C,n)+lam/2*hprime_a2
        A.append(B[n]/(n+1))
        U.append(V[n]/(n+1))
        B.append(rhsA/(n+1))
        V.append(rhsU/(n+1))
        C.append(-sum((A[k]*C[n+1-k] for k in range(1,n+2)),RF(0)))
        T.append(-15*sum((U[k]*T[n+1-k] for k in range(1,n+2)),RF(0))/(15*u0-16))
    theta = x34_cm11_theta(CC)
    power = CC(1)
    b = []
    for an in A:
        b.append(out_CC(an*power))
        power *= theta
    return dict(b=b,theta=out_CC(theta),working_bits=work_bits,
                prec_bits=prec_bits,Nper=Nper,guard_bits=int(guard_bits))


def build_x34_cm11_setup(side_data_path, Nper=1000, prec_bits=1000,
                         guard_bits=None):
    data = _heegner_read_side_data(side_data_path)
    if (QQ(data['quaternion_a']),QQ(data['quaternion_b'])) != (34,-3):
        raise ValueError('Expected quaternion presentation (34,-3).')
    if QQ(data['cm_discriminant']) != -11:
        raise ValueError('Expected CM center discriminant -11.')
    if list(map(QQ,data['cm_element'])) != [QQ(1)/2,QQ(1)/2,-QQ(7)/2,QQ(1)/2]:
        raise ValueError('The supplied recurrence is for the specified -11 lift.')
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out=x34_cm11_series(Nper,prec_bits,guard_bits)
    out.update(tau0=x34_cm11_tau(CC),
               S_matrices=[x34_quaternion_matrix(c,RF) for c in data['side_generators']],
               side_data_path=str(side_data_path),side_data=data,E_label='34a1')
    return out


def audit_x34_series_modularity(setup,Nper=None):
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




# genus_two\x38_cm11_setup.sage
"""Two guarded CM expansions on X(38,1), at D=-11, from Table 2/5.

Sage's original (-19,2) quaternion coordinates are used for CM and AL
exports; the supplied side JSON uses (2,-19), i_new=j_old,
j_new=i_old, k_new=-k_old. Both map to exactly the same real matrices.
"""
from sage.all import *
from pathlib import Path
import json


def x38_quaternion_matrix(coords,RF):
    x,y,z,t=map(QQ,coords)
    s,r=RF(2).sqrt(),RF(19).sqrt()
    return matrix(RF,[[RF(x)+RF(z)*s,(RF(y)-RF(t)*s)*r],
                      [(-RF(y)-RF(t)*s)*r,RF(x)-RF(z)*s]])


def x38_cm11_tau(CC):
    RF=RealField(CC.precision())
    return CC(-2*RF(38).sqrt()/19,RF(209).sqrt()/19)


def x38_cm11_theta(component,CC):
    RF=RealField(CC.precision())
    product_log=sum(ZZ(kronecker(-11,j))*(RF(j)/11).log_gamma() for j in range(1,11))
    pi_omega2=product_log.exp()/22
    coefficient=-RF(9)/19 if component=='a' else RF(13)/19
    return coefficient*pi_omega2*CC(4,RF(22).sqrt())


def x38_cm11_series(component='a',Nper=1600,prec_bits=1000,guard_bits=None):
    assert component in ('a','b')
    Nper,prec_bits=int(Nper),int(prec_bits)
    guard_bits=4*Nper+128 if guard_bits is None else int(guard_bits)
    work_bits=prec_bits+guard_bits
    RF,CC=RealField(work_bits),ComplexField(work_bits)
    d=9 if component=='a' else 13
    lam=-RF(11)/(4*d*d)
    S_coeffs=(RF(33)/2,RF(495)/8,RF(627)/4,RF(627)/8) if component=='a' else (RF(0),RF(44),RF(143),RF(165)/2)
    Su,S0,Sj,Sj2=[v/(d*d) for v in S_coeffs]
    h0,h1,h2=(-328,-472,-192) if component=='a' else (-236,-656,-228)
    A,p,W,U,V,J=[RF(1)],[RF(1)],[RF(1)],[RF(1)],[RF(44)/d],[RF(1)]
    J2,U2,S,HP=[],[],[],[]
    def conv(a,b,n):
        return sum((a[k]*b[n-k] for k in range(n+1)),RF(0))
    for n in range(Nper-1):
        J2.append(conv(J,J,n));U2.append(conv(U,U,n))
        S.append(Su*U[n]+Sj*J[n]+Sj2*J2[n]+(S0 if n==0 else 0))
        HP.append(h1*U[n]+h2*U2[n]+(h0 if n==0 else 0))
        an=conv(p,A,n)/(n+1)
        pn=(conv(p,p,n)/2-conv(S,W,n))/(n+1)
        wn=2*conv(p,W,n)/(n+1)
        un=V[n]/(n+1)
        vn=(conv(p,V,n)+lam/2*conv(W,HP,n))/(n+1)
        jn=-conv(V,J2,n)/(n+1)
        A.append(an);p.append(pn);W.append(wn);U.append(un);V.append(vn);J.append(jn)
    theta=x38_cm11_theta(component,CC)
    outCC=ComplexField(prec_bits);power=CC(1);coeffs=[]
    for an in A:
        coeffs.append(outCC(an*power));power*=theta
    return dict(b=coeffs,theta=outCC(theta),Nper=Nper,prec_bits=prec_bits,
                component=component,working_bits=work_bits,guard_bits=guard_bits)


def build_x38_cm11_setup(component='a',Nper=1600,prec_bits=1000,guard_bits=None):
    side=side_data(38)
    assert (QQ(side['quaternion_a']),QQ(side['quaternion_b']))==(2,-19)
    assert list(map(QQ,side['cm_element']))==[QQ(1)/2,-1,-QQ(1)/2,0]
    out=x38_cm11_series(component,Nper,prec_bits,guard_bits)
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    matrices=[x38_quaternion_matrix([c[0],c[2],c[1],-QQ(c[3])],RF) for c in side['side_generators']]
    out.update(tau0=x38_cm11_tau(CC),S_matrices=matrices,side_data=side,
               side_data_path=None,E_label='38'+component+'1')
    return out


def audit_x38_series_modularity(setup,Nper=None):
    coeffs,tau=setup['b'],setup['tau0'];CC=tau.parent();RF=RealField(CC.precision())
    nper=len(coeffs) if Nper is None else min(int(Nper),len(coeffs));bar=tau.conjugate()
    def mob(g,z):return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])
    def disk(z):return (z-tau)/(z-bar)
    def val(w):
        s=CC(0)
        for b in reversed(coeffs[:nper]):s=s*w+b
        return (1-w)**2*s
    rows=[]
    for i,g in enumerate(setup['S_matrices']):
        end=disk(mob(g.inverse(),tau));w=end/(1+(1-abs(end)**2).sqrt())
        z=(tau-bar*w)/(1-w);gw=disk(mob(g,z))
        left=val(gw);right=(g[1,0]*z+g[1,1])**2/g.det()*val(w)
        error=abs(left-right)/max(RF(1),abs(left),abs(right))
        rows.append(dict(side_index=i,radius=max(abs(w),abs(gw)),residual_bits=CC.precision() if error==0 else -error.log2()))
    return rows




# six_curves\x46_cm163_setup.sage
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
    data = _heegner_read_side_data(side_data_path)
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


# genus_two\x58_cm11_setup.sage
"""X(58,1), both eigenforms at the printed discriminant -11 CM center.

Uses (i^2,j^2)=(2,-29), not the different (58,-3) gamma-data basis.
The recurrence uses guarded ordinary Taylor coefficients.
"""
from sage.all import *
from pathlib import Path
import json,time


def x58_quaternion_matrix(coords,RF):
    x,y,z,t=map(QQ,coords)
    sa,sb=RF(2).sqrt(),RF(29).sqrt()
    return matrix(RF,[[x+y*sa,(z+t*sa)*sb],[-(z-t*sa)*sb,x-y*sa]])


def x58_cm11_tau(CC):
    RF=RealField(CC.precision())
    return CC(-3*RF(58).sqrt()/29,RF(319).sqrt()/29)


def x58_cm11_theta(CC,component):
    RF=RealField(CC.precision())
    lp=sum(ZZ(kronecker(-11,j))*(RF(j)/11).log_gamma() for j in range(1,11))
    omega2=lp.exp()/(2*RF.pi()*11)
    scale=-RF(5)/29 if component=='a' else RF(3)/29
    return scale*RF.pi()*omega2*CC(6,RF(22).sqrt())


def x58_cm11_series(component='a',Nper=1500,prec_bits=1000,guard_bits=None):
    component=str(component)[0]; assert component in ('a','b')
    Nper,prec_bits=int(Nper),int(prec_bits)
    if guard_bits is None: guard_bits=4*Nper+128
    RF=RealField(prec_bits+int(guard_bits));CC=ComplexField(RF.precision())
    if component=='a':
        lam=RF(121)/400;u0=-RF(11);u1=-RF(176)/5
        r0,r1,r2,t1,t2=map(RF,[5,-QQ(144358)/243,-2523,-QQ(229376)/243,QQ(65536)/9])
        aa,bb=RF(1),RF(27)
        hp0,hp1,hp2=map(RF,[-3448,-624,-24])
    else:
        lam=-RF(11)/144;u0=-RF(1)/11;u1=-RF(16)/33
        r0,r1,r2,t1,t2=map(RF,[-QQ(39625)/243,-34,-3,-QQ(163840)/243,-QQ(65536)/243])
        aa,bb=RF(27),RF(1)
        hp0,hp1,hp2=map(RF,[-312,-6896,-20184])
    A,P,W,U,V=[RF(1)],[RF(1)],[RF(1)],[u0],[u1]
    S,T=[1/u0],[1/(aa*u0+bb)]
    S2,T2,Rpart,U2,Hprime=([] for _ in range(5))
    def conv(a,b,n):return sum((a[k]*b[n-k] for k in range(n+1)),RF(0))
    start=time.time()
    for n in range(Nper-1):
        S2.append(conv(S,S,n));T2.append(conv(T,T,n));U2.append(conv(U,U,n))
        Rpart.append((r0 if n==0 else 0)+r1*S[n]+r2*S2[n]+t1*T[n]+t2*T2[n])
        Hprime.append((hp0 if n==0 else 0)+hp1*U[n]+hp2*U2[n])
        an=conv(P,A,n)/(n+1)
        pn=(conv(P,P,n)/2-lam*conv(Rpart,W,n))/(n+1)
        wn=2*conv(P,W,n)/(n+1)
        un=V[n]/(n+1)
        vn=(conv(P,V,n)+lam/2*conv(Hprime,W,n))/(n+1)
        sn=-conv(V,S2,n)/(n+1)
        tn=-aa*conv(V,T2,n)/(n+1)
        A.append(an);P.append(pn);W.append(wn);U.append(un);V.append(vn);S.append(sn);T.append(tn)
    theta=x58_cm11_theta(CC,component)
    power=CC(1);out=ComplexField(prec_bits);coeff=[]
    for an in A:coeff.append(out(an*power));power*=theta
    return dict(b=coeff,theta=out(theta),Nper=Nper,prec_bits=prec_bits,
        guard_bits=int(guard_bits),working_bits=RF.precision(),component=component,
        known_cn=[RealField(prec_bits)(A[n]*factorial(n)) for n in range(min(10,Nper))],
        elapsed=time.time()-start)


def x58_exact_side_data():
    return side_data(58)


def build_x58_cm11_setup(component='a',Nper=1500,prec_bits=1000,guard_bits=None):
    geom=x58_exact_side_data()
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out=x58_cm11_series(component,Nper,prec_bits,guard_bits)
    out.update(tau0=x58_cm11_tau(CC),S_matrices=[x58_quaternion_matrix(c,RF) for c in geom['side_generators']],
               side_data=geom,E_label='58'+component+'1')
    return out


def audit_x58_series_modularity(setup,Nper=None):
    b,tau=setup['b'],setup['tau0'];CC=tau.parent();RF=RealField(CC.precision())
    Nper=min(len(b),len(b) if Nper is None else int(Nper));tbar=tau.conjugate()
    def mob(g,z):return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])
    def disk(z):return (z-tau)/(z-tbar)
    def evaluate(w):
        v=CC(0)
        for an in reversed(b[:Nper]):v=v*w+an
        return (1-w)^2*v
    rows=[]
    for index,g in enumerate(setup['S_matrices']):
        endpoint=disk(mob(g.inverse(),tau));w=endpoint/(1+(1-abs(endpoint)^2).sqrt())
        z=(tau-tbar*w)/(1-w);gz=mob(g,z);gw=disk(gz)
        left=evaluate(gw);right=(g[1,0]*z+g[1,1])^2/g.det()*evaluate(w)
        error=abs(left-right)/max(RF(1),abs(left),abs(right))
        rows.append(dict(side_index=index,radius=max(abs(w),abs(gw)),relative_error=error,
                         residual_bits=CC.precision() if error==0 else -error.log2()))
    return rows
