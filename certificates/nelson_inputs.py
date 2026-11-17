"""Independent exact input checks; run with: sage -python certificates/nelson_inputs.py.

This checks the quaternion, order, CM, and elliptic-curve inputs without
evaluating any Nelson lattice sums.  Requires SageMath and NumPy.
"""
import json
from pathlib import Path

from sage.all import QQ, QQbar, QuaternionAlgebra, QuadraticField, EllipticCurve, matrix
import nelson_cm_values as nelson

ROOT = Path(__file__).resolve().parents[1]


def require(condition, description):
    if not condition:
        raise AssertionError(description)


def algebraic_matrix(a, b, quaternion):
    """Independent algebraic-number version of the stated real embedding."""
    x,y,z,w = [QQ(t) for t in quaternion]
    if a > 0:
        p,q = QQbar(a).sqrt(), QQbar(-b).sqrt()
        return matrix(QQbar, [[x+y*p,(z+w*p)*q],[-(z-w*p)*q,x-y*p]])
    p,q = QQbar(b).sqrt(), QQbar(-a).sqrt()
    return matrix(QQbar, [[x+z*p,(y-w*p)*q],[-(y+w*p)*q,x-z*p]])


def run_certificate(output_dir=None):
    levels = []
    forms = []
    for N, data in nelson.DATA.items():
        algebra = QuaternionAlgebra(QQ,data['a'],data['b'])
        require(algebra.discriminant()==N,f"quaternion discriminant {N}")
        E = nelson.CMPoint(N)
        require(nelson.reduced_discriminant(E.A,E.R)==N,f"order discriminant {N}")
        require(nelson.contains(E.R,(1,0,0,0)),f"order identity {N}")
        require(all(nelson.contains(E.R,E.A.mul(a,b)) for a in E.R for b in E.R),f"order closure {N}")
        require(all((2*a[0]).denominator==1 and E.A.nrd(a).denominator==1 for a in E.R),f"order integral traces/norms {N}")
        dual = nelson.dual(E.A,E.R)
        require(all(E.A.pair(a,b)==int(i==j) for i,a in enumerate(E.R) for j,b in enumerate(dual)),f"trace-dual basis {N}")
        require(nelson.contains(E.R,E.mu),f"CM inclusion {N}")
        require(2*E.mu[0]==1 and 1-4*E.A.nrd(E.mu)==E.D,f"CM discriminant {N}")
        field = QuadraticField(E.D)
        require(field.discriminant()==E.D and field.class_number()==1,f"maximal CM order and class number {N}")
        t = sum(QQ(c)*QQbar(m).sqrt() for m,c in E.tau_exact.d.items())
        M = algebraic_matrix(E.a,E.b,E.mu)
        require(M[1,0]*t*t+(M[1,1]-M[0,0])*t-M[0,1]==0,f"algebraic fixed point {N}")
        require(t.imag()>0,f"exact upper-half-plane choice {N}")
        require(E.tau_exact==nelson.TABLE5_TAU[N],f"tau table equality {N}")
        standard_basis = [tuple(int(i==j) for i in range(4)) for j in range(4)]
        for a in standard_basis:
            for b in standard_basis:
                require(algebraic_matrix(E.a,E.b,E.A.mul(a,b))==algebraic_matrix(E.a,E.b,a)*algebraic_matrix(E.a,E.b,b),f"real matrix embedding multiplication {N}")
        levels.append(dict(N=N,status="PASS",quaternion_discriminant=N,
            order_reduced_discriminant=N,maximal_order=True,order_closed=True,
            trace_dual_verified=True,CM_discriminant=E.D,CM_class_number=1,
            tau=str(E.tau_exact),tau_fixed_point_exact=True,tau_imaginary_part_positive_exact=True,
            real_matrix_embedding_verified=True))
    for form, ainvariants in nelson.CURVES.items():
        N = int(form[:-1])
        curve = EllipticCurve(QQ,ainvariants)
        require(curve.conductor()==N,f"elliptic conductor {form}")
        bad = [int(p) for p in QQ(N).numerator().prime_factors()]
        coefficients = nelson.an_table(ainvariants,100,bad)
        require(coefficients==list(curve.anlist(100)),f"Hecke coefficient check {form}")
        forms.append(dict(form=form,status="PASS",conductor=N,
            Hecke_coefficients_compared=100,independent_reference="Sage EllipticCurve.anlist"))
    report = dict(certificate="Exact input validation for Nelson computations",status="PASS",
        arithmetic="exact QQ and QQbar; Sage quaternion and elliptic-curve arithmetic",
        levels=levels,forms=forms,numerical_lattice_sums_evaluated=False)
    if output_dir is not None:
        path = Path(output_dir)/"nelson_inputs.json"
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    return report


if __name__=="__main__":
    report=run_certificate(ROOT/"results")
    print("PASS: 9 quaternion algebras, maximal orders, exact CM fixed points and upper-half-plane choices.")
    print("PASS: all 12 elliptic conductors and first 100 Hecke coefficients independently checked.")
