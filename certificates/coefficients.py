"""Exact, dependency-free CM expansion certificate for all twelve eigenforms.

Run: python certificates/coefficients.py --count 20 --verify-count 32
API: first_coefficients('14a', 20) returns [c_0,...,c_19], as Fractions.
The theorem is conditional on the input data and algebraic recognition in the
draft. This module certifies consequences of those data; it does not prove
the numerical recognition of the CM period Theta or initial coefficients.
"""
from __future__ import annotations

import argparse
import ast
import json
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
DATA = json.loads((ROOT / "data" / "recurrence_data.json").read_text(encoding="utf-8"))
FORMS = DATA["forms"]


def pad(a, n):
    return list(a[:n]) + [Q(0)] * max(0, n - len(a))


def add(a, b, n):
    return [x + y for x, y in zip(pad(a, n), pad(b, n))]


def scale(a, k, n):
    return [k * x for x in pad(a, n)]


def mul(a, b, n):
    a, b = pad(a, n), pad(b, n)
    return [sum((a[j] * b[i-j] for j in range(i+1)), Q(0)) for i in range(n)]


def inverse(a, n):
    a = pad(a, n)
    if not a[0]:
        raise ZeroDivisionError("The denominator vanishes at the expansion center")
    b = [1 / a[0]]
    for i in range(1, n):
        b.append(-sum((a[j] * b[i-j] for j in range(1, i+1)), Q(0)) / a[0])
    return b


def power(a, exponent, n):
    if exponent < 0:
        return power(inverse(a, n), -exponent, n)
    result = pad([Q(1)], n)
    while exponent:
        if exponent & 1:
            result = mul(result, a, n)
        exponent >>= 1
        if exponent:
            a = mul(a, a, n)
    return result


def derivative(a):
    return [j * a[j] for j in range(1, len(a))]


def polynomial(coefficients, u, n):
    value = [Q(0)] * n
    for c in reversed(coefficients):
        value = mul(value, u, n)
        value[0] += Q(c)
    return value


@lru_cache(maxsize=None)
def expression_tree(expression):
    return ast.parse(expression, mode="eval").body


def rational_series(expression, u, n):
    """Evaluate a whitelisted rational expression exactly in Q[[q]]/(q^n)."""
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return pad([Q(node.value)], n)
        if isinstance(node, ast.Name) and node.id == "u":
            return pad(u, n)
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return scale(visit(node.operand), Q(-1), n)
            if isinstance(node.op, ast.UAdd):
                return visit(node.operand)
        if isinstance(node, ast.BinOp):
            left = visit(node.left)
            if isinstance(node.op, ast.Pow):
                exponent = ast.literal_eval(node.right)
                if type(exponent) is not int:
                    raise ValueError("Only integer exponents are allowed")
                return power(left, exponent, n)
            right = visit(node.right)
            if isinstance(node.op, ast.Add):
                return add(left, right, n)
            if isinstance(node.op, ast.Sub):
                return add(left, scale(right, Q(-1), n), n)
            if isinstance(node.op, ast.Mult):
                return mul(left, right, n)
            if isinstance(node.op, ast.Div):
                return mul(left, inverse(right, n), n)
        raise ValueError("Unsupported expression node: " + ast.dump(node))
    return visit(expression_tree(expression))


def _validate_request(form, count):
    if form not in FORMS:
        raise ValueError(f"Unknown form {form!r}; choose from {', '.join(FORMS)}")
    if type(count) is not int or count < 0:
        raise ValueError("count must be a nonnegative integer")


def ordinary_series(form, count):
    """Return ordinary coefficient lists for A(q),u(q), through q^(count-1).

    For step k, the right-hand sides use only coefficients through k+1.
    Thus dividing the known coefficient by (k+1)(k+2) uniquely appends the
    next coefficients, over Q, with no numerical root or branch choices.
    """
    _validate_request(form, count)
    d = FORMS[form]
    A = [Q(1), Q(1)]
    u = [Q(d["u0"]), Q(d["u1"])]
    L = Q(d["L"])
    Hp = derivative([Q(x) for x in d["H"]])
    for k in range(max(0, count - 2)):
        n = k + 1
        invA = inverse(A, n)
        Ap, up = derivative(A), derivative(u)
        A2 = mul(A, A, n)
        rhsA = add(scale(mul(mul(Ap, Ap, n), invA, n), Q(3, 2), n),
                   scale(mul(rational_series(d["R"], u, n), mul(A2, A, n), n), -L, n), n)
        rhsu = add(mul(mul(Ap, invA, n), up, n),
                   scale(mul(A2, polynomial(Hp, u, n), n), L / 2, n), n)
        denominator = Q((k+1) * (k+2))
        A.append(rhsA[k] / denominator)
        u.append(rhsu[k] / denominator)
    return A[:count], u[:count]


def first_coefficients(form, count):
    """Return exactly count normalized coefficients c_0,...,c_(count-1)."""
    A, _ = ordinary_series(form, count)
    return [factorial(n) * a for n, a in enumerate(A)]


def appendix_expansion(form, count):
    """Return the appendix normalization and exact c_n and c_n/n! arrays."""
    A, _ = ordinary_series(form, count)
    d = FORMS[form]
    return {"form": form, "D": d["D"], "tau": d["tau"],
            "theta_over_pi_Omega_D_squared": d["theta_over_pi_Omega_D_squared"],
            "formula": "f(z)=(1-w)^2*sum(c_n*(Theta*w)^n/n!); w=(z-tau)/(z-conjugate(tau))",
            "c": [factorial(n) * a for n, a in enumerate(A)],
            "ordinary_A": A}


def w_coefficients(form, count, theta=None):
    """Coefficients of f(z), with (1-w)^2 multiplied out.

    With theta=None, return polynomials in the formal symbol Theta as maps
    {degree: Fraction}; otherwise substitute the caller's value of Theta.
    The n-th coefficient is a_n*Theta^n-2*a_(n-1)*Theta^(n-1)
    +a_(n-2)*Theta^(n-2), where a_n=c_n/n! and negative indices mean zero.
    """
    A, _ = ordinary_series(form, count)
    out = []
    for n in range(count):
        p = {n: A[n]}
        if n >= 1:
            p[n-1] = -2 * A[n-1]
        if n >= 2:
            p[n-2] = A[n-2]
        out.append(p if theta is None else sum(a * theta**k for k, a in p.items()))
    return out


def verify_form(form, count=32):
    """Check all printed data and independent undivided residual identities."""
    if count < 10:
        raise ValueError("Verification requires count >= 10")
    A, u = ordinary_series(form, count)
    d = FORMS[form]
    c = [factorial(k) * a for k, a in enumerate(A)]
    for key in ("assumed_c", "appendix_c"):
        expected = [Q(x) for x in d[key]]
        assert c[:len(expected)] == expected, (form, key, c[:len(expected)], expected)
    # Check original equations, not the divided update formulas.
    # count known terms imply second-derivative residual through q^(count-3).
    n = count - 2
    L = Q(d["L"])
    Ap, up = derivative(A), derivative(u)
    App, upp = derivative(Ap), derivative(up)
    R = rational_series(d["R"], u, n)
    H = polynomial(d["H"], u, n)
    Hp = polynomial(derivative([Q(x) for x in d["H"]]), u, n)
    r1 = add(add(mul(A, App, n), scale(mul(Ap, Ap, n), Q(-3, 2), n), n),
             scale(mul(R, power(A, 4, n), n), L, n), n)
    r2 = add(add(mul(A, upp, n), scale(mul(Ap, up, n), Q(-1), n), n),
             scale(mul(power(A, 3, n), Hp, n), -L/2, n), n)
    # This first integral is additional to the second-order recursion.
    integral = add(mul(up, up, n), scale(mul(mul(A, A, n), H, n), -L, n), n)
    for label, residual in (("Schwarzian", r1), ("coordinate", r2), ("first integral", integral)):
        assert all(x == 0 for x in residual), (form, label, residual)
    # Exercise the original-w API with an exact nontrivial substitution and
    # compare direct multiplication by (1-w)^2.
    theta = Q(7, 3)
    expected_w = mul([Q(1), Q(-2), Q(1)], [a * theta**k for k, a in enumerate(A)], count)
    assert w_coefficients(form, count, theta) == expected_w
    return {"form": form, "status": "PASS", "coefficients_computed": count,
            "assumed_coefficients_matched": len(d["assumed_c"]),
            "appendix_coefficients_matched": len(d["appendix_c"]),
            "undivided_equation_residual_coefficients_zero": n,
            "first_integral_residual_coefficients_zero": n,
            "original_w_prefactor_check": "PASS"}


def run(count=20, verify_count=32, output_dir=None):
    if count < 1:
        raise ValueError("Display count must be positive")
    verify_count = max(verify_count, count, 10)
    output_dir = Path(output_dir) if output_dir else ROOT / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    results, expansions = [], {}
    for form in FORMS:
        results.append(verify_form(form, verify_count))
        expansions[form] = appendix_expansion(form, count)
        print(f"{form}: PASS; c_0,...,c_{count-1} = " + ", ".join(map(str, expansions[form]["c"])), flush=True)
    report = {"status": "PASS", "scope": "exact conditional recurrence consequences",
              "normalization": DATA["normalization"], "forms": results,
              "elapsed_seconds": round(perf_counter()-start, 3),
              "expansions": expansions}
    (output_dir / "coefficients.json").write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    lines = ["# First normalized coefficients of all twelve eigenforms", "",
             f"Each list contains exactly {count} entries: c_0 through c_{count-1}.", "",
             "The displayed coefficients occur in f(z)=(1-w)^2 sum_n c_n (Theta w)^n/n!, with w=(z-tau)/(z-conjugate(tau)).",
             "They are exact rational consequences of Table 2. The recognition of the input coefficients and Theta is a separate hypothesis.", ""]
    for form, item in expansions.items():
        lines.extend([f"## {form}", "", f"tau = `{item['tau']}`", "",
                      f"Theta/(pi Omega_D^2) = `{item['theta_over_pi_Omega_D_squared']}`", "",
                      "| n | c_n |", "|---:|:---|"])
        lines.extend(f"| {n} | {value} |" for n, value in enumerate(item["c"]))
        lines.append("")
    (output_dir / "coefficients_first20.md").write_text("\n".join(lines), encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=20)
    parser.add_argument("--verify-count", type=int, default=32)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    run(args.count, args.verify_count, args.output_dir)
