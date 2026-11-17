"""Exact cross-module consistency certificate for the twelve draft forms.

Requires a fresh recovered_parameters.json export from the executed Sage
notebook and the Appendix C certificate output. Missing dependencies fail;
no cached manuscript value is substituted for recovered parameters.
No Nelson numerical sum or under-200 image computation is invoked.
"""
from __future__ import annotations

import importlib.util
import json
from hashlib import sha256
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
U = s.Symbol("u")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact(value):
    if hasattr(value, "numerator") and hasattr(value, "denominator"):
        return s.Rational(value.numerator, value.denominator)
    return s.sympify(value, locals={"u": U, "sqrt": s.sqrt})


def equal(left, right, message):
    difference = s.simplify(s.cancel(s.expand(exact(left)-exact(right))))
    require(difference == 0, message + ": " + str(difference))


def polynomial(coefficients):
    return sum(exact(c)*U**n for n, c in enumerate(coefficients))


def rad_expression(rad):
    return sum(exact(coefficient)*s.sqrt(radicand) for radicand, coefficient in rad.d.items())


def normalize_form(label):
    label = str(label)
    if label.endswith("1") and label[-2] in "ab":
        label = label[:-1]
    return label


def normalize_seed_list(values):
    if isinstance(values, dict):
        values = {int(n): v for n, v in values.items()}
        require(set(values) == set(range(max(values)+1)), "Notebook coefficients must be contiguous from c_0")
        return [values[n] for n in range(max(values)+1)]
    require(isinstance(values, list), "Notebook coefficients must be a list or indexed mapping")
    return values


def run():
    required = [
        "data/recurrence_data.json",
        "results/recovered_parameters.json",
        "results/appendix_c.json",
        "certificates/families.py",
        "certificates/coefficients.py",
        "certificates/nelson_cm_values.py",
    ]
    missing = [name for name in required if not (ROOT/name).is_file()]
    require(not missing, "Missing prerequisite certificate output/source: " + ", ".join(missing))
    recurrence = json.loads((ROOT/required[0]).read_text(encoding="utf-8"))["forms"]
    recovered_raw = json.loads((ROOT/required[1]).read_text(encoding="utf-8"))
    recovered = {normalize_form(label): value for label, value in recovered_raw.items()}
    require(len(recovered) == len(recovered_raw), "Duplicate notebook form labels after normalization")
    appendix = json.loads((ROOT/required[2]).read_text(encoding="utf-8"))
    require(appendix["status"] == "PASS", "Appendix C prerequisite did not pass")
    families = load_module("certificate_integration_families", ROOT/required[3])
    coefficients = load_module("certificate_integration_coefficients", ROOT/required[4])
    nelson = load_module("certificate_integration_nelson", ROOT/required[5])
    forms = set(recurrence)
    require(len(forms) == 12, "Expected exactly twelve recurrence forms")
    for name, dataset in (("recovered notebook parameters", recovered),
                          ("family fibers", families.FIBERS),
                          ("Nelson input coefficients", nelson.TABLE4),
                          ("Nelson Theta", nelson.TABLE5_THETA),
                          ("coefficient engine", coefficients.FORMS)):
        require(set(dataset) == forms, name + " does not cover the same twelve forms")
    appendix_rows = {(int(r["N"]), int(r["D"])): r for r in appendix["coordinate_rows"]}
    require(len(appendix_rows) == len(appendix["coordinate_rows"]), "Duplicate Appendix C coordinate rows")
    output_rows = []
    for form, d in recurrence.items():
        level = int(form[:-1])
        recovered_form = recovered[form]
        require(set(("H", "R", "L", "u0", "u1", "c")).issubset(recovered_form), "Missing notebook fields for " + form)
        H, R = polynomial(d["H"]), exact(d["R"])
        H_derived = families.derive_H(form)[0]
        equal(H, H_derived, form + " recurrence H versus model-derived H")
        family = families.derive_family(form)
        fixed = exact(family["fixed_part"])
        basis = [exact(b) for b in family["accessory_basis"]]
        parameters = s.symbols("accessory0:" + str(len(basis)))
        residual_numerator = s.fraction(s.cancel(R-fixed-sum(a*b for a,b in zip(parameters,basis))))[0]
        equations = s.Poly(residual_numerator, U).all_coeffs()
        solutions = list(s.linsolve(equations, parameters))
        require(len(solutions) == 1 and all(value.is_Rational for value in solutions[0]), form + " R is not a unique rational member of the derived family")
        equal(R, fixed+sum(a*b for a,b in zip(solutions[0],basis)), form + " family membership residual")
        recovered_H = polynomial(recovered_form["H"]) if isinstance(recovered_form["H"], list) else exact(recovered_form["H"])
        equal(H, recovered_H, form + " H versus recovered notebook")
        equal(R, recovered_form["R"], form + " R versus recovered notebook")
        for name in ("L", "u0", "u1"):
            equal(d[name], recovered_form[name], form + " " + name + " versus recovered notebook")
        seeds = [exact(c) for c in d["assumed_c"]]
        notebook_seeds = [exact(c) for c in normalize_seed_list(recovered_form["c"])]
        require(len(notebook_seeds) >= len(seeds), form + " notebook lacks assumed seed coefficients")
        require(notebook_seeds[:len(seeds)] == seeds, form + " assumed seeds versus notebook")
        generated = list(map(exact, coefficients.first_coefficients(form, len(notebook_seeds))))
        require(notebook_seeds == generated, form + " notebook input coefficients versus recurrence output")
        nelson_seeds = nelson.TABLE4[form]
        require(set(nelson_seeds) == set(range(2, len(seeds))), form + " Nelson seed index coverage")
        require(all(exact(nelson_seeds[n]) == seeds[n] for n in nelson_seeds), form + " Nelson seeds versus recurrence")
        require(nelson.DATA[level]["D"] == d["D"] and form in nelson.DATA[level]["forms"], form + " Nelson discriminant/form mapping")
        equal(rad_expression(nelson.TABLE5_TAU[level]), d["tau"], form + " Nelson tau versus appendix metadata")
        equal(rad_expression(nelson.TABLE5_THETA[form]), d["theta_over_pi_Omega_D_squared"], form + " Nelson Theta algebraic factor versus appendix metadata")
        coordinate_row = appendix_rows.get((level, d["D"]))
        require(coordinate_row is not None and coordinate_row["status"] == "PASS", form + " missing passed Appendix C center")
        coordinate_key = "u" if "u" in coordinate_row else "u_"+form[-1]
        equal(coordinate_row[coordinate_key], d["u0"], form + " Appendix C center versus recurrence")
        equal(exact(d["u1"])**2, exact(d["L"])*H.subs(U, exact(d["u0"])), form + " first-integral initial condition")
        denominator = s.denom(s.cancel(R))
        require(denominator.subs(U, exact(d["u0"])) != 0, form + " R has a pole at recurrence center")
        output_rows.append({"form": form, "status": "PASS",
            "model_H_matches_recurrence_and_notebook": True,
            "R_matches_recovered_notebook_and_derived_family": True,
            "accessory_basis": list(map(str, basis)), "accessory_values": list(map(str, solutions[0])),
            "L_u0_u1_match_recovered_notebook": True,
            "notebook_coefficients_checked": len(notebook_seeds),
            "Nelson_assumed_coefficients_checked": len(nelson_seeds),
            "Nelson_tau_exact_identity": True, "Nelson_Theta_factor_exact_identity": True,
            "Appendix_C_center_exact_identity": True, "first_integral_initial_condition": True})
        print(form + ": PASS cross-module identities", flush=True)
    report = {"status": "PASS", "form_count": len(output_rows),
              "scope": "Exact consistency of independent certificate inputs and freshly recovered notebook parameters; no new numerical evaluation.",
              "prerequisite_order": ["execute notebooks and export recovered_parameters.json", "run Appendix C", "run this integration check"],
              "dependencies": [{"path": name, "sha256": sha256((ROOT/name).read_bytes()).hexdigest()} for name in required],
              "conditionality": "Agreement of exact data and formal consequences does not prove analytic identification of numerically recognized initial coefficients or Theta.",
              "forms": output_rows}
    (ROOT/"results"/"integration.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    try:
        run()
    except Exception as error:
        output = ROOT/"results"/"integration.json"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps({"status": "FAIL", "error": str(error)}, indent=2)+"\n", encoding="utf-8")
        raise
