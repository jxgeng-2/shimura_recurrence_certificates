"""Exact arithmetic certificate for every Appendix C coordinate.

Run: python certificates/appendix_c.py
Requires SymPy. Published arithmetic/geometric theorems remain explicit inputs;
this script certifies the coordinate algebra and weighted moduli identities.
"""
import json
from math import gcd, isqrt
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
j, x, T, z = s.symbols("j x T z")
Q = s.Rational
WEIGHTS = (2, 3, 5, 6)


def require(value, message):
    if not value:
        raise AssertionError(message)


def eq(a, b):
    return s.cancel(a-b) == 0


def at(expr, value, variable=j):
    if value == s.oo:
        value = s.limit(expr, variable, s.oo)
    else:
        value = s.cancel(expr).subs(variable, value)
    return s.oo if value in (s.zoo, -s.oo) else value


def reduced_forms(D):
    """Primitive reduced positive definite forms of negative discriminant D."""
    out = []
    for a in range(1, isqrt((-D)//3)+1):
        for b in range(-a, a+1):
            if (b*b-D) % (4*a):
                continue
            c = (b*b-D)//(4*a)
            if a > c or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            out.append((a, b, c))
    return out


def quadratic_symbol(D, p):
    if D % p == 0:
        return 0
    if p == 2:
        return 1 if D % 8 in (1,7) else -1
    return 1 if pow(int(D)%int(p),(int(p)-1)//2,int(p)) == 1 else -1


def parametrizations():
    """Lin--Yang Table 4, expanded arXiv 1807.00466, pp.44--46.

    s_k are the coordinates of equation (15), of weighted degrees 2,3,5,6.
    This retains complete polynomials, so substitution is independently rerun.
    """
    return {
        21: (
            (5*j**2-16*j-5)**2,
            7*(j-1)*(j+1)*(27*j**4-144*j**3+394*j**2-624*j+603),
            2**12*(j-1)*(j+1)*(j+3)**2*(j**2-5)**2,
            2**12*(j**12+6*j**11-9*j**10-108*j**9-28*j**8+708*j**7
                    +406*j**6-2016*j**5+17*j**4+3846*j**3-2637*j**2-2436*j+6346)),
        33: (
            288*j**6-360*j**5-295*j**4+240*j**3+238*j**2+120*j+25,
            (j-1)*(j+1)*(324*j**8-4536*j**7+9216*j**6+4284*j**5+3373*j**4
                          +4704*j**3+3350*j**2+924*j+121),
            2**6*(j-1)*j**2*(j+1)*(3*j-1)**2*(3*j**2+1)**2*(3*j**3+7*j**2-3*j+1)**2,
            2**4*(6561*j**20+26244*j**19+275562*j**18+26244*j**17+279693*j**16
                  +863136*j**15-190440*j**14-779904*j**13+700018*j**12-123384*j**11
                  -72484*j**10+20328*j**9+35570*j**8-28416*j**7+12888*j**6
                  -3936*j**5+1133*j**4-300*j**3+74*j**2-12*j+1)),
        46: (
            (j-1)**2*(j+1)**2*(6*j**3+5*j**2-9*j+5)**2,
            (j-1)**2*(j+1)**2*(324*j**10-504*j**9+226*j**8+252*j**7-7*j**6
                              +1260*j**5-2214*j**4-252*j**3+1681*j**2-756*j+98),
            3**5*(j-1)**4*j**2*(j+1)**4*(2*j+1)**2*(j**2-j-1)**2
                 *(2*j**2-1)**2*(3*j**2-1)**2,
            3**5*(j-1)**4*(j+1)**4*(192*j**22-192*j**21+1152*j**20+2304*j**19
                 -12*j**18-7164*j**17-13633*j**16+8806*j**15+33381*j**14-2310*j**13
                 -38998*j**12-6064*j**11+26610*j**10+8116*j**9-10904*j**8-4788*j**7
                 +2510*j**6+1532*j**5-242*j**4-258*j**3-11*j**2+18*j+3))}


def coordinate_changes():
    # P_N are the exact draft models. F_N are branch polynomials from
    # Lin--Yang Table 7; proportionality certifies the marked base change.
    models = {
        21: -7*x**4+94*x**2-343,
        33: -3*x**4-10*x**2-243,
        34: -3*x**4+26*x**3-53*x**2-26*x-3,
        46: -x**4+45*x**2-512}
    branches = {
        21: 7*j**4-94*j**2+343,
        33: 243*j**4+10*j**2+3,
        34: 27*j**4+136*j**3+122*j**2-136*j+27,
        46: 64*j**4-45*j**2+8}
    changes = {21:7/j, 33:1/j, 34:(5+3*j)/(3-5*j), 46:8*j}
    y_scales = {21:7/j**2, 33:1/j**2, 34:17/(3-5*j)**2, 46:s.Integer(8)}
    expected_ratios = {21:-49/j**4, 33:-1/j**4, 34:-289/(3-5*j)**4, 46:s.Integer(-64)}
    results = []
    for N in changes:
        ratio = s.factor(models[N].subs(x, changes[N])/branches[N])
        require(eq(ratio, expected_ratios[N]), f"N={N}: branch transformation")
        require(eq(ratio,-y_scales[N]**2), f"N={N}: explicit y-coordinate lift")
        numerator, denominator = s.fraction(changes[N])
        require(s.degree(numerator,j) <= 1 and s.degree(denominator,j) <= 1,
                f"N={N}: expected Mobius change")
        require(s.diff(changes[N],j) != 0, f"N={N}: nonconstant change")
        results.append({"N":N, "x(j)":str(changes[N]), "P(x(j))/F(j)":str(ratio),
                        "y/Y for auxiliary presentation Y^2=-F(j)":str(y_scales[N]),
                        "status":"PASS exact branch polynomial identity"})
    # Elkies equations (68), (76): full model isomorphisms over Q.
    P14 = -x**4+13*x**2-128
    E14 = -16*z**4+13*z**2-8
    require(eq(P14.subs(x,4*z),16*E14), "Elkies 14 model x=4s, y=4v")
    E15 = -(3*z**2+1)*(z**2+27)
    P15 = -3*x**4-82*x**2-27
    require(eq(P15.subs(x,z),E15), "Elkies 15 model x=s, y=v")
    results[0:0] = [
        {"N":14,"full_model_change":"x=4s; y=4v; T=-s^2; u=x^2=-16T", "status":"PASS"},
        {"N":15,"full_model_change":"x=s; y=v; T=-3s^2; u=x^2=-T/3", "status":"PASS"}]
    return results


def coordinate_rows():
    u34 = s.cancel((5+3*j)/(3-5*j)-(3-5*j)/(5+3*j))
    transformations = {14:-16*j,15:-j/3,21:49/j**2,33:1/j**2,34:u34,46:64*j**2}
    for N in (21,33,34,46):
        involution=-1/j if N==34 else -j
        require(eq(transformations[N].subs(j,involution),transformations[N]),
                f"N={N}: Atkin-Lehner quotient invariance")
    # Source values are independent inputs, not inferred from target u values.
    rows = [
        (14,-11,-1,16,"Elkies section 5.1"),
        (14,-4,s.oo,s.oo,"Elkies section 5.1"),
        (15,-7,-27,9,"Elkies Table 6"),
        (15,-3,s.oo,s.oo,"Elkies Table 6"),
        (21,-43,Q(5,3),Q(441,25),"Lin-Yang/Errthum intersection certificate"),
        (21,-4,1,49,"Lin-Yang Table 6"),
        (33,-67,Q(5,3),Q(9,25),"Lin-Yang/Errthum intersection certificate"),
        (33,-3,s.oo,0,"Lin-Yang Table 6"),
        (33,-4,1,1,"Lin-Yang Table 6"),
        (34,-11,1,Q(-15,4),"Lin-Yang Table 6"),
        (34,-3,s.oo,Q(16,15),"Lin-Yang Table 6"),
        (46,-163,Q(7,9),Q(3136,81),"Lin-Yang/Errthum intersection certificate"),
        (46,-4,s.oo,s.oo,"Lin-Yang Table 6"),
        (46,-3,1,64,"Lin-Yang Table 6")]
    result = []
    for N,D,v,target,source in rows:
        actual = at(transformations[N],v)
        require(actual == target,f"N={N}, D={D}: coordinate value")
        result.append({"N":N,"D":D,"published_parameter":str(v),"u":str(actual),"source":source,"status":"PASS"})
    # Guo--Yang coordinates: s=x^2 for 26/38; s=-x^2/19 for 58.
    genus2 = [
        (26,-11,1,(1,1)), (26,-52,0,(0,s.oo)), (26,-8,s.oo,(s.oo,0)),
        (38,-11,1,(1,1)), (38,-19,0,(0,s.oo)), (38,-4,s.oo,(s.oo,0)),
        (58,-11,Q(11,19),(-11,Q(-1,11))), (58,-232,0,(0,s.oo)),
        (58,-8,s.oo,(s.oo,0)), (58,-3,Q(27,19),(-27,Q(-1,27)))]
    for N,D,v,target in genus2:
        ua = j if N != 58 else -19*j
        actual = (at(ua,v),at(1/ua,v))
        require(actual == target,f"N={N}, D={D}: reciprocal coordinates")
        result.append({"N":N,"D":D,"published_s":str(v),"u_a":str(actual[0]),"u_b":str(actual[1]),
                       "source":{26:"Guo-Yang Example 33",38:"Guo-Yang Table 11",58:"Guo-Yang Table 16"}[N],"status":"PASS"})
    require(at(transformations[21],s.oo)==0,"marked CM_21(-7) value")
    # The other two marked reciprocal choices occur among the 24 rows above.
    require(len(result)==24,"Appendix C coverage")
    return result


def intersections():
    polys=parametrizations()
    centers=[(21,-43,Q(5,3),10,21),(33,-67,Q(5,3),20,231),(46,-163,Q(7,9),1870,27531)]
    # Errthum Table 3, two independently transcribed factored norms.
    norms={
        -43:(Q(3**7*7**4,2**10*5**6),Q(19**2*37**2*43,2**10*5**6)),
        -67:(Q(3**7*7**4*11**4,2**16*5**6),Q(13**2*43**2*61**2*67,2**16*5**6)),
        -163:(Q(3**11*7**4*19**4*23**4,2**10*5**6*11**6*17**6),
              Q(13**2*67**2*109**2*139**2*157**2*163,2**10*5**6*11**6*17**6))}
    result=[]
    for N,D,j0,a,b in centers:
        raw=tuple(s.cancel(p.subs(j,j0)) for p in polys[N])
        target=tuple(map(s.Integer,(a*a,b*b,a*a*b*b,a**6+b**4)))
        lam=s.cancel((raw[1]/target[1])/(raw[0]/target[0]))
        for w,r,t in zip(WEIGHTS,raw,target):
            require(eq(r,lam**w*t),f"N={N}, weighted s_{w} substitution")
        q=Q(b**4,a**6)
        S6=(q,q*q,q**3,q**3*(q+1))
        mu=Q(a**4,b*b)
        for w,t,v in zip(WEIGHTS,target,S6):
            require(eq(t,mu**w*v),f"N={N}: S6 equivalence degree {w}")
        norm_t,norm_one_minus_t=norms[D]
        require(norm_t == 27*q/16,f"D={D}: first Errthum norm")
        require(norm_one_minus_t == 1+27*q/16,f"D={D}: second Errthum norm")
        # If t=v+iw, |t|^2=A^2 and |1-t|^2=(1+A)^2 force v=-A,w=0.
        re_t=s.cancel((1+norm_t**2-norm_one_minus_t**2)/2)
        im_t_squared=s.cancel(norm_t**2-re_t**2)
        require(re_t==-norm_t and im_t_squared==0,f"D={D}: unique sign")
        require(eq(-16*re_t/27,q),f"D={D}: J=-16 t6/27")
        forms=reduced_forms(D)
        require(len(forms)==1,f"D={D}: class number one")
        require(s.ntheory.primetest.is_square(abs(D)) is False,"negative fundamental discriminant")
        require(s.isprime(-D) and D%4==1,"maximal order discriminant")
        local={int(p):quadratic_symbol(D,int(p)) for p in s.factorint(N)}
        require(all(v!=1 for v in local.values()),f"D={D}: quaternion local embedding condition")
        require(raw[2]!=0, "nonzero s5 / nondecomposable surface")
        result.append({"N":N,"D":D,"j":str(j0),"a":a,"b":b,"q":str(q),
                       "S_N_raw":[str(v) for v in raw],"raw_equals_lambda_power_times_target":str(lam),
                       "target":[str(v) for v in target],"weighted_scale_S6_to_target":str(mu),
                       "Errthum_abs_t6":str(norm_t),"Errthum_abs_1_minus_t6":str(norm_one_minus_t),
                       "t6":str(re_t),"class_number":1,"reduced_forms":forms,
                       "local_kronecker_symbols":local,"status":"PASS exact arithmetic; CM conclusion uses moduli/commutant argument"})
    return result


def run():
    changes=coordinate_changes()
    rows=coordinate_rows()
    common=intersections()
    h={str(D):{"h":len(reduced_forms(D)),"reduced_forms":reduced_forms(D)}
       for D in sorted({r["D"] for r in rows})}
    require(h["-52"]["h"]==h["-232"]["h"]==2,"zero-fiber class numbers")
    return {"certificate":"Appendix C","status":"PASS","coordinate_row_count":len(rows),
            "coordinate_changes":changes,"coordinate_rows":rows,"three_center_intersections":common,
            "class_numbers":h,"scope":"Exact arithmetic over Q. Literature CM identifications and geometric moduli/commutant argument are explicit theorem inputs; no numerical recognition used.",
            "source_records":"data/appendix_c_sources.json","proof_notes":"docs/appendix_c.md"}


if __name__=="__main__":
    result=run()
    path=ROOT/"results"/"appendix_c.json"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("PASS Appendix C: 24 table rows, 6 coordinate changes, 3 four-coordinate weighted identities, 3 two-norm sign determinations.")
    print(path)
