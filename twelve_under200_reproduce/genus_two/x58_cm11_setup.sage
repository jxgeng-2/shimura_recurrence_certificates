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
    return _heegner_read_side_data('x58_cm11_side_data.json')


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
