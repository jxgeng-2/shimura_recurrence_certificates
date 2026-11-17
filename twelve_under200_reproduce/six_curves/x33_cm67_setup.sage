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

def build_x33_cm67_setup(side_data_path='x33_cm67_side_data.json',Nper=1000,prec_bits=1000,guard_bits=None):
    with open(side_data_path,encoding='utf-8-sig') as f: data=json.load(f)
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
