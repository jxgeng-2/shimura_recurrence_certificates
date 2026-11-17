"""Symbolic certificate for the identities behind the paper's recurrences.

This checks identities with arbitrary symbolic jets, not chosen numerical
examples. It does not turn the paper's geometric theorems into a formal proof.
Run: python certificates/general_identities.py
"""
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
Q = s.Rational


def check_zero(expression, name, results):
    residual = s.cancel(expression)
    if residual != 0:
        raise AssertionError((name, residual))
    results.append({"identity": name, "residual": "0", "status": "PASS"})


def schwarzian(f, q):
    return s.diff(f,q,3)/s.diff(f,q)-Q(3,2)*(s.diff(f,q,2)/s.diff(f,q))**2


def run(output_dir=None):
    results = []
    q, z = s.symbols("q z")
    g0,g1,g2,g3,h1,h2,h3 = s.symbols("g0 g1 g2 g3 h1 h2 h3")
    g = g0+g1*z+g2*z**2/2+g3*z**3/6
    h = h1*q+h2*q**2/2+h3*q**3/6
    # Arbitrary third jets determine a Schwarzian. h(0)=0 loses no generality
    # because the expansion center in g's coordinate may be translated.
    direct = schwarzian(g.subs(z,h),q).subs(q,0)
    chain = schwarzian(g,z).subs(z,0)*h1**2+schwarzian(h,q).subs(q,0)
    check_zero(direct-chain,"Schwarzian chain rule for arbitrary third jets",results)

    a,b,c,d = s.symbols("a b c d")
    mobius = (a*q+b)/(c*q+d)
    check_zero(schwarzian(mobius,q),"Schwarzian of an arbitrary nondegenerate Mobius map",results)
    tau,taubar,theta = s.symbols("tau taubar theta")
    inverse_disk = (tau-taubar*q/theta)/(1-q/theta)
    check_zero(schwarzian(inverse_disk,q),"Schwarzian of z as a function of normalized disk q",results)
    w=s.symbols("w")
    inverse_w = (tau-taubar*w)/(1-w)
    check_zero((1-w)**2*s.diff(inverse_w,w)-(tau-taubar),
               "Cancellation of modular-form prefactor against dz/dw",results)

    A,Ap,App,U,Up,Upp,L,R,H,Hp,Hpp,lam = s.symbols("A Ap App U Up Upp L R H Hp Hpp lambda")
    # t'=lambda*A, t''=lambda*A', t'''=lambda*A''.
    t_schwarzian = lam*App/(lam*A)-Q(3,2)*(lam*Ap/(lam*A))**2
    chain_with_zero_left = R*lam**2*A**2+t_schwarzian
    recurrence = A*App-Q(3,2)*Ap**2+L*R*A**4
    check_zero((A**2*chain_with_zero_left).subs(lam**2,L)-recurrence,
               "dt=lambda*A dq gives the undivided Schwarzian recurrence",results)

    # Chain differentiation from u_t^2=H(u); lambda^2=L.
    ut = s.symbols("u_t",nonzero=True)
    coordinate_second_derivative = Hp*lam**2*A**2/2+ut*lam*Ap
    check_zero(coordinate_second_derivative.subs(ut,Up/(lam*A)).subs(lam**2,L)
               -(Ap*Up/A+L*A**2*Hp/2),
               "Coordinate equation from u_t^2=H and t_q=lambda*A",results)

    # d/dq[(u'/A)^2-L H(u)] = (2u'/A^3)*(A u''-A' u'-L A^3 H'/2).
    first_integral_derivative=2*Up*Upp/A**2-2*Up**2*Ap/A**3-L*Hp*Up
    coordinate_residual=A*Upp-Ap*Up-L*A**3*Hp/2
    check_zero(first_integral_derivative-2*Up*coordinate_residual/A**3,
               "Exact propagation of the first integral (u'/A)^2-LH(u)",results)
    F=Up**2-L*A**2*H
    Fprime=2*Up*Upp-2*L*A*Ap*H-L*A**2*Hp*Up
    check_zero(Fprime.subs(Upp,Ap*Up/A+L*A**2*Hp/2)-2*Ap*F/A,
               "Equivalent error evolution F'=2(A'/A)F",results)

    # The conversion formula used for the three published Schwarzians.
    utt=Hp/2
    uttt=Hpp*ut/2
    schwarzian_u_t=uttt/ut-Q(3,2)*(utt/ut)**2
    check_zero(schwarzian_u_t.subs(ut**2,H)-(Hpp/2-Q(3,8)*Hp**2/H),
               "Quotient correction {u,t}=H''/2-3(H')^2/(8H)",results)
    phi1,phi2,phi3,Qphi=s.symbols("phi1 phi2 phi3 Qphi")
    phi_schwarzian=phi3/phi1-Q(3,2)*(phi2/phi1)**2
    converted=H*(2*Qphi*phi1**2+phi_schwarzian)+schwarzian_u_t.subs(ut**2,H)
    check_zero(converted-(2*H*Qphi*phi1**2+H*phi_schwarzian+Hpp/2-Q(3,8)*Hp**2/H),
               "Published Q(T) to R(u) conversion with arbitrary coordinate jets",results)

    # A nonelliptic uniformization contributes an arbitrary holomorphic
    # Schwarzian S0+S1*v+... . Its contribution starts below both tested poles.
    v,m,aa,bb,cc,S0,S1=s.symbols("v m a0 b0 c0 S0 S1")
    k=aa+bb*v+cc*v**2
    logh=m/v+s.diff(k,v)/k
    scht=s.diff(logh,v)-logh**2/2
    scaled_R=s.cancel(((S0+S1*v)-scht)*v**2/k**2)
    leading=s.cancel(scaled_R.subs(v,0))
    subleading=s.cancel(s.diff(scaled_R,v).subs(v,0))
    check_zero(leading-m*(m+2)/(2*aa**2),
               "Generic nonelliptic local leading polar coefficient",results)
    check_zero(subleading+m*(m+1)*bb/aa**3,
               "Generic nonelliptic local second polar coefficient",results)
    e=s.symbols("e",positive=True)
    scaled_elliptic=s.cancel(((1-1/e**2)/(2*v**2)+S0/v-scht)*v**2/k**2)
    check_zero(scaled_elliptic.subs(v,0)-((m+1)**2-1/e**2)/(2*aa**2),
               "Generic elliptic local leading polar coefficient",results)

    # Coefficient-reconstruction jets in Section 4, with c_n=A^(n)(0).
    c2,c3,c4,c5=s.symbols("c2 c3 c4 c5")
    Aq=1+q+c2*q**2/2+c3*q**3/6+c4*q**4/24+c5*q**5/120
    target=(Q(3,2)*s.diff(Aq,q)**2-Aq*s.diff(Aq,q,2))/Aq**4
    check_zero(target.subs(q,0)-(Q(3,2)-c2),"Reconstruction jet LR(u0)=3/2-c2",results)
    check_zero(s.diff(target,q).subs(q,0)-(6*c2-c3-6),
               "Reconstruction jet LR'(u0)u1=6c2-c3-6",results)

    # The recursive updates uniquely solve the two undivided equations.
    # Independent arbitrary jets verify their next-unknown coefficient:
    # coefficient q^n of A*A'' is (n+1)(n+2)*a_(n+2) plus known jets.
    triangular=[]
    for n in range(7):
        ac=s.symbols("a0:"+str(n+3))
        uc=s.symbols("u0:"+str(n+3))
        polyA=1+sum(ac[k]*q**k for k in range(1,n+3))
        polyu=sum(uc[k]*q**k for k in range(n+3))
        rA=s.Poly(s.expand(polyA*s.diff(polyA,q,2)-Q(3,2)*s.diff(polyA,q)**2),q).nth(n)
        ru=s.Poly(s.expand(polyA*s.diff(polyu,q,2)-s.diff(polyA,q)*s.diff(polyu,q)),q).nth(n)
        denominator=(n+1)*(n+2)
        check_zero(s.diff(rA,ac[n+2])-denominator,f"Triangular A coefficient at q^{n}",results)
        check_zero(s.diff(ru,uc[n+2])-denominator,f"Triangular u coefficient at q^{n}",results)
        check_zero(s.diff(ru,ac[n+2]),f"Coordinate equation has no new A coefficient at q^{n}",results)
        triangular.append({"q_degree":n,"new_coefficient_multiplier":denominator})

    report={"certificate":"General Schwarzian and recurrence identities","status":"PASS",
            "arithmetic":"exact rational identities with symbolic jets",
            "identity_checks":results,"checks_count":len(results),"triangular_checks":triangular,
            "nonvanishing_conditions":["Mobius determinant nonzero", "local derivatives in a Schwarzian nonzero", "A(0)=1", "local leading differential coefficient nonzero"],
            "scope":"Algebraic identities in Sections 2-4. Geometry, uniformization and analytic convergence are theorem inputs. The displayed finite triangular checks illustrate the general coefficient argument; they are not asserted as a finite proof for every index."}
    output_dir=Path(output_dir) if output_dir else ROOT/"results"
    output_dir.mkdir(parents=True,exist_ok=True)
    (output_dir/"general_identities.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    return report


if __name__=="__main__":
    report=run()
    print(f"PASS: {report['checks_count']} exact general-identity and symbolic coefficient checks.")
