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


def build_x14_cm11_setup(side_data_path='x14_cm11_side_data.json',Nper=1400,prec_bits=1000,guard_bits=None):
    with open(side_data_path,encoding='utf-8-sig') as f:
        data=json.load(f)
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
