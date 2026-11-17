"""Exact, small sanity certificate for the FOUR Section 5 Table 4 claims.

No CM period integration, high-order recurrence, algdep, or under-200 sweep
is run. Polynomial properties and matching exact quaternion inputs are
checked. Analytic CM-image identification is deliberately not claimed.
Uses only the Python standard library.
"""
from __future__ import annotations

import ast
import json
import re
from fractions import Fraction as Q
from functools import reduce
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def primes_up_to(bound):
    return [n for n in range(2, bound+1) if all(n % d for d in range(2, isqrt(n)+1))]


def prime_factors(n):
    return [p for p in primes_up_to(n) if n % p == 0]


def trim(a):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    return a


def polynomial_mod(a, b, p):
    a, b = trim([x % p for x in a]), trim([x % p for x in b])
    if not b:
        raise ZeroDivisionError("zero polynomial")
    inverse_lead = pow(b[-1], -1, p)
    while len(a) >= len(b):
        c, offset = a[-1] * inverse_lead % p, len(a)-len(b)
        for j in range(len(b)):
            a[j+offset] = (a[j+offset]-c*b[j]) % p
        a = trim(a)
    return a


def polynomial_gcd(a, b, p):
    while b:
        a, b = b, polynomial_mod(a, b, p)
    if not a:
        return []
    c = pow(a[-1], -1, p)
    return [c*x % p for x in a]


def polynomial_power_mod(a, exponent, f, p):
    def product(a, b):
        result = [0] * (len(a)+len(b)-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                result[i+j] = (result[i+j]+x*y) % p
        return polynomial_mod(result, f, p)
    result = [1]
    while exponent:
        if exponent & 1:
            result = product(result, a)
        exponent >>= 1
        if exponent:
            a = product(a, a)
    return result


def subtract_x(a, p):
    a = list(a) + [0] * max(0, 2-len(a))
    a[1] = (a[1]-1) % p
    return trim(a)


def irreducibility_certificate(coeffs):
    """Find/check a finite-field irreducibility witness (Rabin criterion).

    Primitive f in Z[x] is irreducible over Q if a degree-preserving
    reduction is irreducible. For degree n, verify x^(p^n)=x mod f and
    gcd(f,x^(p^(n/r))-x)=1 for every prime divisor r of n.
    """
    n = len(coeffs)-1
    if n == 1:
        return {"method": "nonconstant linear polynomial", "degree": 1}
    for p in primes_up_to(1000):
        if coeffs[-1] % p == 0:
            continue
        f = [x % p for x in coeffs]
        if subtract_x(polynomial_power_mod([0, 1], p**n, f, p), p):
            continue
        checks = []
        for r in prime_factors(n):
            rem = subtract_x(polynomial_power_mod([0, 1], p**(n//r), f, p), p)
            g = polynomial_gcd(f, rem, p)
            if g != [1]:
                break
            checks.append({"prime_divisor_of_degree": r, "gcd": [1]})
        else:
            return {"method": "degree-preserving irreducible reduction over F_p",
                    "prime": p, "coefficients_mod_p_ascending": f,
                    "frobenius_p_power_degree_remainder": [0], "gcd_checks": checks}
    raise AssertionError("No small-prime irreducibility witness found")


def reduced_forms(D):
    """Enumerate primitive reduced positive binary forms, including boundaries."""
    result = []
    for a in range(1, isqrt((-D)//3)+1):
        for b in range(-a, a+1):
            if (b*b-D) % (4*a):
                continue
            c = (b*b-D)//(4*a)
            if c < a or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            result.append((a, b, c))
    return result


def squarefree(n):
    return all(n % (d*d) != 0 for d in range(2, isqrt(n)+1))


def fundamental(D):
    return (D % 4 == 1 and squarefree(abs(D))) or (D % 4 == 0 and (D//4) % 4 in (2, 3) and squarefree(abs(D//4)))


def kronecker_at_prime(D, p):
    if D % p == 0:
        return 0
    if p == 2:
        return 1 if D % 8 in (1, 7) else -1
    return 1 if pow(D % p, (p-1)//2, p) == 1 else -1


def rational_literal(node):
    """Read only exact literal data, never execute supplied Sage source."""
    if isinstance(node, ast.Constant) and type(node.value) in (int, str):
        return Q(node.value) if type(node.value) is int else node.value
    if isinstance(node, (ast.List, ast.Tuple)):
        return [rational_literal(x) for x in node.elts]
    if isinstance(node, ast.Dict):
        return {rational_literal(k): rational_literal(v) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -rational_literal(node.operand)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        return rational_literal(node.left) / rational_literal(node.right)
    raise ValueError("Nonliteral source input: " + ast.dump(node))


def exact_source_data(path):
    names = {"CM_DB", "CM_ALGEBRA_I2", "CM_ALGEBRA_J2", "CM_ORDER_BASIS",
             "CM_AL_COORDS", "CM_CLASS_RECORDS", "CM_ALL_RECORDS", "CM_SELECTED_DS", "CM_RAMIFIED_PRIMES"}
    result = {}
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in names:
                result[name] = rational_literal(node.value)
    assert set(result) == names
    return result


def quaternion_product(x, y, a, b):
    x0, x1, x2, x3 = x
    y0, y1, y2, y3 = y
    return [x0*y0+a*x1*y1+b*x2*y2-a*b*x3*y3,
            x0*y1+x1*y0-b*x2*y3+b*x3*y2,
            x0*y2+x2*y0+a*x1*y3-a*x3*y1,
            x0*y3+x3*y0+x1*y2-x2*y1]


def quaternion_norm(v, a, b):
    return v[0]**2-a*v[1]**2-b*v[2]**2+a*b*v[3]**2


def basis_coordinates(v, basis):
    rows = [[basis[col][row] for col in range(4)] + [v[row]] for row in range(4)]
    for k in range(4):
        pivot = next(i for i in range(k, 4) if rows[i][k])
        rows[k], rows[pivot] = rows[pivot], rows[k]
        c = rows[k][k]
        rows[k] = [x/c for x in rows[k]]
        for i in range(4):
            if i != k:
                c = rows[i][k]
                rows[i] = [x-c*y for x, y in zip(rows[i], rows[k])]
    return [row[-1] for row in rows]


def check_quaternion_inputs(claim):
    N, D, h, operator = (claim[k] for k in ("N", "D", "h", "operator"))
    path = ROOT / "twelve_under200_reproduce" / "cm_data" / f"x{N}_cm_data.sage"
    d = exact_source_data(path)
    assert d["CM_DB"] == N and D in d["CM_SELECTED_DS"] and operator in d["CM_AL_COORDS"]
    a, b, basis = d["CM_ALGEBRA_I2"], d["CM_ALGEBRA_J2"], d["CM_ORDER_BASIS"]
    ps = list(map(int, d["CM_RAMIFIED_PRIMES"]))
    assert ps == prime_factors(N)
    chi = [kronecker_at_prime(D, p) for p in ps]
    assert 1 not in chi
    orientations = reduce(lambda x, y: x*y, (1-c for c in chi), 1)
    classes = [r for r in d["CM_CLASS_RECORDS"] if r[0] == D]
    all_records = [r for r in d["CM_ALL_RECORDS"] if r[0] == D]
    assert len(classes) == h and len(all_records) == h*orientations
    assert {(int(r[6]), int(r[7])) for r in all_records} == {(o, c) for o in range(1, orientations+1) for c in range(1, h+1)}
    for r in classes+all_records:
        assert int(r[1]) == h and list(map(int, r[3])) == chi
        mu = r[4]
        trace, norm = 2*mu[0], quaternion_norm(mu, a, b)
        assert trace*trace-4*norm == D
        assert trace.denominator == norm.denominator == 1
        assert all(x.denominator == 1 for x in basis_coordinates(mu, basis))
        # A fundamental discriminant and integral generator certify the full
        # maximal quadratic order Z[mu] at these odd discriminants.
        assert trace == 1 and norm == Q(1-D, 4)
    w = d["CM_AL_COORDS"][operator]
    norm = quaternion_norm(w, a, b)
    assert norm == int(operator[1:]) > 0
    winv = [w[0]/norm]+[-x/norm for x in w[1:]]
    normalizer_columns = []
    for v in basis:
        conjugate = quaternion_product(quaternion_product(w, v, a, b), winv, a, b)
        coords = basis_coordinates(conjugate, basis)
        assert all(x.denominator == 1 for x in coords)
        normalizer_columns.append([int(x) for x in coords])
    # Check the inverse inclusion too, not just containment in the order.
    for v in basis:
        conjugate = quaternion_product(quaternion_product(winv, v, a, b), w, a, b)
        assert all(x.denominator == 1 for x in basis_coordinates(conjugate, basis))
    return {"source": str(path.relative_to(ROOT)), "source_sha256": sha256(path.read_bytes()).hexdigest(),
            "chi_at_ramified_primes": chi, "orientation_count": orientations,
            "class_representatives_checked": len(classes), "all_orientation_representatives_checked": len(all_records),
            "mu_trace": "1", "mu_norm": str(Q(1-D, 4)), "mu_discriminant": D,
            "mu_integral_in_supplied_order": True, "operator_coordinates": list(map(str, w)),
            "operator_norm": str(norm), "operator_normalizes_supplied_order_both_directions": True,
            "operator_conjugation_columns_in_order_basis": normalizer_columns}


def run(output_dir=None):
    payload = json.loads((ROOT / "data" / "cm_claims.json").read_text(encoding="utf-8"))
    source_root = ROOT / "twelve_under200_reproduce"
    sources = [p for p in source_root.rglob("*") if p.suffix in (".sage", ".m", ".json", ".txt")]
    source_texts = {str(p.relative_to(ROOT)): p.read_text(encoding="utf-8") for p in sources}
    fingerprint = re.search(r"EXPECTED_TABLE_SHA256\s*=\s*['\"]([0-9a-f]{64})", (source_root / "export_table.sage").read_text(encoding="utf-8")).group(1)
    rows = []
    for claim in payload["claims"]:
        coeffs = list(map(int, claim["coefficients"]))
        assert coeffs[-1] > 0 and reduce(gcd, coeffs) == 1
        assert fundamental(claim["D"]) and 3 <= abs(claim["D"]) < 200
        forms = reduced_forms(claim["D"])
        assert len(forms) == claim["h"] == len(coeffs)-1
        irreducibility = irreducibility_certificate(coeffs)
        inputs = check_quaternion_inputs(claim)
        occurrences = [name for name, content in source_texts.items() if str(coeffs[-1]) in content]
        assert not occurrences, ("Source bundle changed; review polynomial provenance", occurrences)
        row = dict(claim, status="PASS_EXACT_SANITY", primitive=True,
                   reduced_positive_quadratic_forms=forms, irreducibility=irreducibility,
                   quaternion_input_checks=inputs, polynomial_source="manuscript Table 4",
                   leading_coefficient_occurrences_in_regeneration_source=occurrences,
                   analytic_CM_image_identification="NOT_RERUN", exact_H_D_point_lifting="NOT_RUN")
        rows.append(row)
        print(f"{claim['E']}, D={claim['D']}, {claim['operator']}: PASS exact sanity; analytic identification NOT RERUN", flush=True)
    report = {"status": "PASS_EXACT_SANITY", "full_under200_run": False,
              "analytic_CM_image_identification": "NOT_RERUN",
              "scope": "Four manuscript polynomials: primitive irreducibility, degree/class number, admissibility, exact selected CM quaternion records and operator normalizers.",
              "source_bundle_whole_table_fingerprint": fingerprint,
              "fingerprint_recomputed": False,
              "fingerprint_limit": "The fingerprint cannot recover or independently check selected polynomial rows without regenerating the whole table; it is recorded only as source provenance.",
              "rows": rows}
    output_dir = Path(output_dir) if output_dir else ROOT / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "cm_claims.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    run()
