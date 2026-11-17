"""Exact model-to-accessory-family certificate for Table 3 (all 12 forms).

This derives H from (model, differential, quotient coordinate), derives the
prescribed polar coefficients from local Schwarzian calculations, and derives
the accessory basis from the descended divisor.  The elliptic labels and their
locations are geometric inputs certified separately in Appendix C.

Run: python3 certificates/families.py [--output FILE]
Requires SymPy.  No floating point computation or external service is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as s

x, u, v = s.symbols("x u v")
C, B, B0, B1 = s.symbols("C B B0 B1")
Q = s.Rational

# Model data are transcribed independently of the claimed H and R families.
MODELS = {
    14: -x**4 + 13*x**2 - 128,
    15: -3*x**4 - 82*x**2 - 27,
    21: -7*x**4 + 94*x**2 - 343,
    33: -3*x**4 - 10*x**2 - 243,
    34: -3*x**4 + 26*x**3 - 53*x**2 - 26*x - 3,
    46: -x**4 + 45*x**2 - 512,
    26: -2*x**6 + 19*x**4 - 24*x**2 - 169,
    38: -16*x**6 - 59*x**4 - 82*x**2 - 19,
    58: -2*x**6 - 78*x**4 - 862*x**2 - 1682,
}

# a=None denotes infinity; each tuple is (a, e, r, m).
FIBERS = {
    "14a": [(None, 2, 2, 0)],
    "15a": [(None, 3, 2, 0)],
    "21a": [(Q(49), 2, 1, 0)],
    "33a": [(Q(0), 3, 2, 0), (Q(1), 2, 1, 0)],
    "34a": [(Q(16, 15), 3, 1, 0)],
    "46a": [(None, 2, 2, 0), (Q(64), 3, 1, 0)],
    "26a": [(Q(0), 1, 2, 1)],
    "26b": [(Q(0), 1, 2, 1)],
    "38a": [(Q(0), 1, 2, 1), (None, 2, 2, 0)],
    "38b": [(Q(0), 2, 2, 1)],
    "58a": [(Q(0), 1, 2, 1), (Q(-27), 3, 1, 0)],
    "58b": [(Q(0), 1, 2, 1), (Q(-1, 27), 3, 1, 0)],
}

# Independent transcriptions of the columns being certified.
EXPECTED_H = {
    "14a": 4*u*(-u**2 + 13*u - 128),
    "15a": 4*u*(-3*u**2 - 82*u - 27),
    "21a": 4*u*(-7*u**2 + 94*u - 343),
    "33a": 4*u*(-3*u**2 - 10*u - 243),
    "34a": (u**2+4)*(-3*u**2 + 26*u - 59),
    "46a": 4*u*(-u**2 + 45*u - 512),
    "26a": 4*(-2*u**3 + 19*u**2 - 24*u - 169),
    "26b": 4*(-169*u**3 - 24*u**2 + 19*u - 2),
    "38a": 4*(-16*u**3 - 59*u**2 - 82*u - 19),
    "38b": 4*(-19*u**3 - 82*u**2 - 59*u - 16),
    "58a": 4*(-2*u**3 - 78*u**2 - 862*u - 1682),
    "58b": 4*(-1682*u**3 - 862*u**2 - 78*u - 2),
}
EXPECTED_FAMILIES = {
    "14a": C-3*u/8,
    "15a": C-4*u/3,
    "21a": C+B/(u-49)-921984/(u-49)**2,
    "33a": C+B/(u-1)-384/(u-1)**2-108/u,
    "34a": C+B/(u-Q(16,15))-Q(1336336,16875)/(u-Q(16,15))**2,
    "46a": C+B/(u-64)-196608/(u-64)**2-3*u/8,
    "26a": C+B/u-Q(507,2)/u**2,
    "26b": C+B/u-3/u**2,
    "38a": C+B/u-Q(57,2)/u**2-6*u,
    "38b": C+B/u-30/u**2,
    "58a": C+B0/u+B1/(u+27)-2523/u**2+Q(65536,9)/(u+27)**2,
    "58b": C+B0/u+B1/(u+Q(1,27))-3/u**2-Q(65536,177147)/(u+Q(1,27))**2,
}
EXPECTED_DIMENSIONS = dict(zip(FIBERS, (1,1,2,2,2,2,2,2,2,2,3,3)))


def check(condition, message):
    """Checks are not disabled by python -O."""
    if not condition:
        raise AssertionError(message)


def zero(expr):
    return s.cancel(expr) == 0


def schwarzian(f, variable):
    return s.diff(f, variable, 3)/s.diff(f, variable) - Q(3,2)*(s.diff(f, variable, 2)/s.diff(f, variable))**2


def kappa(e):
    return (1-Q(1,e**2))/2


def local_schwarzian_certificate():
    """Differentiate z=v^(1/e), t=alpha*v^(m+1)/(m+1).

    Lower analytic terms do not change the displayed highest polar term.
    This verifies the local coefficient used for every fiber type here.
    """
    alpha = s.Symbol("alpha", nonzero=True)
    rows = []
    for m, e in ((0,2), (0,3), (1,1), (1,2)):
        z = v**Q(1,e)
        t = alpha*v**(m+1)/(m+1)
        computed = s.simplify((schwarzian(z,v)-schwarzian(t,v))/s.diff(t,v)**2)
        expected = ((m+1)**2-Q(1,e**2))/(2*alpha**2)*v**(-2*m-2)
        check(zero(computed-expected), f"local Schwarzian m={m},e={e}")
        rows.append({"m":m, "e":e, "leading_expression":str(computed)})
    return rows


def derive_H(form):
    level = int(form[:-1])
    model = MODELS[level]
    genus = (s.degree(model,x)-2)//2
    omega_numerator = x if genus == 2 and form.endswith("a") else s.Integer(1)
    coordinate = x-1/x if level == 34 else (x**-2 if form.endswith("b") else x**2)
    pulled_H = s.cancel(s.diff(coordinate,x)**2*model/omega_numerator**2)
    if level == 34:
        # An exact resultant descent avoids assuming the claimed quartic H.
        z = s.Symbol("z")
        numerator = s.fraction(s.cancel(z-pulled_H))[0]
        relation = x**2-u*x-1
        remainder = s.rem(numerator,relation,x)
        solved = s.solve(s.Poly(remainder,x).all_coeffs(),z)
        H = s.factor(solved[z])
    else:
        polynomial = s.Poly(model,x)
        P = sum(coeff*u**(power[0]//2) for power,coeff in polynomial.terms())
        H = 4*u*P if genus == 1 else (4*P if form.endswith("a") else s.expand(4*u**3*P.subs(u,1/u)))
    check(zero(pulled_H-H.subs(u,coordinate)), f"H pullback {form}")
    check(zero(H-EXPECTED_H[form]), f"H table equality {form}")
    check(s.degree(s.gcd(model,s.diff(model,x)),x)==0, f"smooth model {form}")

    # Verify generators preserve the model, coordinate, and eigenline.
    sigma_x, sigma_y_factor = (-1/x,-1/x**2) if level==34 else (-x,s.Integer(1))
    check(zero(model.subs(x,sigma_x)-sigma_y_factor**2*model), f"sigma model {form}")
    check(zero(coordinate.subs(x,sigma_x)-coordinate), f"sigma u {form}")
    sigma_omega = s.cancel(omega_numerator.subs(x,sigma_x)*s.diff(sigma_x,x)/(sigma_y_factor*omega_numerator))
    check(sigma_omega in (-1,1), f"sigma eigenline {form}")
    # Hyperelliptic involution fixes x,u and sends omega to -omega.
    # [C(x,y):C(x)]=2; deg(x -> u)=2 gives quotient degree 4.
    coordinate_num, coordinate_den = s.fraction(s.cancel(coordinate))
    degree_x_to_u = max(s.degree(coordinate_num,x),s.degree(coordinate_den,x))
    check(degree_x_to_u==2, f"quotient degree {form}")
    return H, coordinate, omega_numerator, model, int(genus), int(sigma_omega)


def derive_family(form):
    H, coordinate, omega_num, model, genus, sigma_sign = derive_H(form)
    fixed = s.Integer(0)
    basis = [s.Integer(1)]
    local_rows = []
    delta = []
    for a,e,r,m in FIBERS[form]:
        delta_order = (2*m+int(e>1))//r
        check(delta_order in (0,1), f"unsupported divisor coefficient {form}")
        row = {"a": "infinity" if a is None else str(a), "e":e,"r":r,"m":m,"delta_order":delta_order}
        if a is None:
            check((r,m)==(2,0), f"infinity type {form}")
            check(zero(coordinate.subs(x,1/v)-v**-2), f"infinity ramification {form}")
            dv_dt_squared = s.cancel(v**4*model.subs(x,1/v)/omega_num.subs(x,1/v)**2)
            leading = s.limit(dv_dt_squared,v,0)
            check(leading != 0 and leading.is_finite, f"infinity regular omega {form}")
            coefficient = kappa(e)*leading
            fixed += coefficient*u
            row.update({"coefficient_of_u":str(coefficient),"dv_dt_squared_at_infinity":str(leading)})
        elif (m,r)==(0,1):
            check(H.subs(u,a)!=0, f"unramified elliptic {form}")
            coefficient = kappa(e)*H.subs(u,a)
            fixed += coefficient/(u-a)**2
            row.update({"double_pole_coefficient":str(coefficient),"H_at_a":str(H.subs(u,a))})
        elif (m,r)==(0,2):
            check(H.subs(u,a)==0 and s.diff(H,u).subs(u,a)!=0, f"ramified elliptic {form}")
            check(a==0 and coordinate==x**2 and model.subs(x,0)!=0, f"local branch model {form}")
            coefficient = kappa(e)*s.diff(H,u).subs(u,a)/4
            fixed += coefficient/(u-a)
            row.update({"simple_pole_coefficient":str(coefficient),"H_prime_at_a":str(s.diff(H,u).subs(u,a))})
        elif (m,r)==(1,2):
            check(a==0 and genus==2, f"zero fiber model {form}")
            # omega^2=(alpha^2*v^2+O(v^4))*dv^2 at the zero.
            if form.endswith("a"):
                local_omega_squared = v**2/model.subs(x,v)
                local_u = coordinate.subs(x,v)
            else:
                local_omega_squared = s.cancel(v**-4/model.subs(x,1/v))
                local_u = coordinate.subs(x,1/v)
            check(zero(local_u-v**2), f"zero ramification {form}")
            alpha_squared = s.limit(local_omega_squared/v**2,v,0)
            check(alpha_squared != 0 and alpha_squared.is_finite, f"simple differential zero {form}")
            check(s.diff(s.cancel(local_omega_squared/v**2),v).subs(v,0)==0, f"zero odd term {form}")
            check(zero(H.subs(u,a)-4/alpha_squared), f"zero H coefficient {form}")
            coefficient = kappa(2*e)*H.subs(u,a)
            fixed += coefficient/(u-a)**2
            row.update({"double_pole_coefficient":str(coefficient),"H_at_a":str(H.subs(u,a)),"alpha_squared":str(alpha_squared)})
        else:
            raise AssertionError(f"Unexpected fiber {form}: {(a,e,r,m)}")
        if delta_order:
            check(a is not None, f"unexpected accessory infinity {form}")
            delta.append(str(a))
            basis.append(1/(u-a))
        local_rows.append(row)

    parameters = [C] + ([B] if len(basis)==2 else ([B0,B1] if len(basis)==3 else []))
    family = fixed+sum(par*b for par,b in zip(parameters,basis))
    check(zero(family-EXPECTED_FAMILIES[form]), f"Table 3 family {form}")
    dimension = 1+len(delta)
    check(dimension==EXPECTED_DIMENSIONS[form], f"Table 3 dimension {form}")
    denominator = s.prod(u-s.sympify(a) for a in delta)
    polynomial_basis = [s.Poly(s.cancel(b*denominator),u) for b in basis]
    coefficient_matrix = s.Matrix([[p.nth(k) for p in polynomial_basis] for k in range(dimension)])
    check(coefficient_matrix.rank()==dimension, f"accessory basis independence {form}")

    # Confirm exactly the requested principal parts and no extra finite poles.
    for row in local_rows:
        if row["a"]=="infinity":
            check(zero(s.limit(family/u,u,s.oo)-s.sympify(row["coefficient_of_u"])), f"infinity principal part {form}")
        else:
            a = s.sympify(row["a"])
            if "double_pole_coefficient" in row:
                check(zero(s.limit((u-a)**2*family,u,a)-s.sympify(row["double_pole_coefficient"])), f"double principal part {form}")
            else:
                check(zero(s.limit((u-a)*family,u,a)-s.sympify(row["simple_pole_coefficient"])), f"simple principal part {form}")
    return {
        "form":form,"status":"PASS","genus":genus,
        "model_y_squared":str(model),"coordinate_u":str(coordinate),
        "differential":f"({omega_num}) dx/y", "sigma_differential_sign":sigma_sign,
        "quotient_degree":4,"H":str(s.factor(H)),
        "pole_fibers":local_rows,"delta_support":delta,
        "accessory_basis":[str(b) for b in basis],"dimension":dimension,
        "fixed_part":str(fixed),"family":str(family),
        "table_H_residual":"0","table_family_residual":"0",
    }


def run_certificate(output_dir=None):
    report = {
        "certificate":"Table 3: model-derived H and all accessory families",
        "arithmetic":"exact rational functions over Q (SymPy)",
        "scope":"Conditional only on the cited model, eigendifferential and elliptic-fiber identifications; no numerical coefficient recognition is used.",
        "local_schwarzian_checks":local_schwarzian_certificate(),
        "forms":[derive_family(form) for form in FIBERS],
        "status":"PASS","forms_checked":len(FIBERS),
    }
    if output_dir is not None:
        path = Path(output_dir)/"families_certificate.json"
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parents[1]/"data"/"families_certificate.json")
    args = parser.parse_args()
    report = run_certificate()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    for row in report["forms"]:
        print(f"{row['form']}: d={row['dimension']}; R(u) = {row['family']}")
    print(f"PASS: all {report['forms_checked']} H functions, principal parts, accessory bases and Table 3 families.")
    print(f"Certificate written to {args.output}")


if __name__=="__main__":
    main()
