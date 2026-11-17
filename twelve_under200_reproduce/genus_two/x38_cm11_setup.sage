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
    with open('x38_cm11_side_data.json',encoding='utf-8-sig') as stream: side=json.load(stream)
    assert (QQ(side['quaternion_a']),QQ(side['quaternion_b']))==(2,-19)
    assert list(map(QQ,side['cm_element']))==[QQ(1)/2,-1,-QQ(1)/2,0]
    out=x38_cm11_series(component,Nper,prec_bits,guard_bits)
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    matrices=[x38_quaternion_matrix([c[0],c[2],c[1],-QQ(c[3])],RF) for c in side['side_generators']]
    out.update(tau0=x38_cm11_tau(CC),S_matrices=matrices,side_data=side,
               side_data_path='x38_cm11_side_data.json',E_label='38'+component+'1')
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


