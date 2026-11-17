"""Guarded ordinary-Taylor recurrences from the supplied Table 2 and Table 5.

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
    side=_heegner_read_side_data('x26_cm11_side_data.json')
    assert (QQ(side['quaternion_a']),QQ(side['quaternion_b']))==(26,-5)
    assert list(map(QQ,side['cm_element']))==[QQ(1)/2,1,-QQ(7)/2,-QQ(1)/2]
    out=x26_cm11_series(form,Nper,prec_bits,guard_bits)
    RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    out.update(tau0=x26_cm11_tau(CC),S_matrices=[x26_quaternion_matrix(c,RF) for c in side['side_generators']],
        gamma_matrices=[x26_quaternion_matrix(c,RF) for c in X26_GAMMA_COORDS],
        side_data=side,side_data_path='x26_cm11_side_data.json',E_label=form+'1')
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
