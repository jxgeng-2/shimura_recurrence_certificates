"""
Generic CM-image and algebraic-recognition package.

Public functions
----------------

make_magma_cm_data(...)
    Package the variables printed by the Magma CM-point script.

generate_cm_image_table(...)
    Compute [w(tau_D)] - [tau_D], recognize its x-coordinate with algdep,
    and write CSV, LaTeX, Sage, and checkpoint files.

recover_recognized_points_over_subfields(...)
    For every recognized row, recover exact points over Q(x,y), record the
    separate fields Q(x) and Q(y), and verify their inclusions in H_D.

Before calling generate_cm_image_table, load the file that defines
heegner_image_from_series, or pass that function explicitly with the optional
heegner_function keyword.
"""

from sage.all import (
    QQ, ZZ, ComplexField, EllipticCurve, Integer, NumberField, PolynomialRing,
    QuadraticField, RealField, latex, load, log, matrix, save,
)
from sage.arith.misc import algdep

import csv
import os
import re
import warnings


CM_IMAGE_TABLE_PACKAGE_VERSION = 1


def make_magma_cm_data(
    cm_points,
    algebra_i2,
    algebra_j2,
    discriminant=None,
    ramified_primes=None,
):
    """
    Package the output of the Magma CM-point script.

    The Magma coordinates must be relative to

        1, i, j, i*j,  with i^2 = algebra_i2 and j^2 = algebra_j2.
    """
    return {
        "cm_points": list(cm_points),
        "algebra_i2": algebra_i2,
        "algebra_j2": algebra_j2,
        "discriminant": discriminant,
        "ramified_primes": (
            None if ramified_primes is None else list(ramified_primes)
        ),
    }


def magma_cm_data_from_namespace(namespace=None):
    """
    Package variables pasted directly from the Magma calculator.

    Example:

        cm_data = magma_cm_data_from_namespace(globals())
    """
    ns = globals() if namespace is None else namespace
    required = ("cm_points", "CM_ALGEBRA_I2", "CM_ALGEBRA_J2")
    missing = [name for name in required if name not in ns]
    if missing:
        raise ValueError(
            "Missing pasted Magma variables: %s" % ", ".join(missing)
        )
    return make_magma_cm_data(
        ns["cm_points"],
        ns["CM_ALGEBRA_I2"],
        ns["CM_ALGEBRA_J2"],
        discriminant=ns.get("CM_DB"),
        ramified_primes=ns.get("CM_RAMIFIED_PRIMES"),
    )


def _cm_first(mapping, names, default=None):
    for name in names:
        if name in mapping:
            return mapping[name]
    return default


def _normalize_cm_data(cm_data):
    if not isinstance(cm_data, dict):
        raise TypeError(
            "cm_data must be a dictionary made by make_magma_cm_data(...)."
        )

    points = _cm_first(cm_data, ("cm_points", "points", "CM_POINTS"))
    if points is None:
        raise ValueError("cm_data has no 'cm_points' list.")

    normalized = []
    for record in points:
        if isinstance(record, dict):
            if "D" not in record or "mu_coords" not in record:
                raise ValueError(
                    "Every CM record needs keys 'D' and 'mu_coords'."
                )
            normalized.append(dict(record))
        elif isinstance(record, (list, tuple)) and len(record) == 2:
            normalized.append({
                "D": record[0],
                "mu_coords": record[1],
            })
        else:
            raise TypeError(
                "A CM record must be a dictionary or [D, mu_coords]."
            )

    grouped = {}
    for record in normalized:
        D = Integer(record["D"])
        record["D"] = D
        grouped.setdefault(D, []).append(record)

    for D, records in grouped.items():
        records.sort(key=lambda rec: Integer(rec.get("orbit_index", 10**9)))
        count = Integer(len(records))
        declared = {
            Integer(rec["h"]) for rec in records if "h" in rec
        }
        if len(declared) > 1:
            raise ValueError("Inconsistent class numbers in the records for D=%s." % D)
        h = next(iter(declared)) if declared else count
        if h != count:
            warnings.warn(
                "D=%s declares h=%s but supplies %s CM records. Recognition "
                "will use degree bound 2h and orbit diagnostics from the "
                "available records." % (D, h, count)
            )
        for index, record in enumerate(records, start=1):
            record.setdefault("h", h)
            record.setdefault("orbit_index", Integer(index))

    return {
        "records": normalized,
        "records_by_D": grouped,
        "algebra_i2": _cm_first(
            cm_data,
            ("algebra_i2", "CM_ALGEBRA_I2", "i2"),
        ),
        "algebra_j2": _cm_first(
            cm_data,
            ("algebra_j2", "CM_ALGEBRA_J2", "j2"),
        ),
        "discriminant": _cm_first(
            cm_data,
            ("discriminant", "CM_DB", "DB"),
        ),
        "ramified_primes": _cm_first(
            cm_data,
            ("ramified_primes", "CM_RAMIFIED_PRIMES"),
        ),
    }


def _cm_exact_rational(value):
    if isinstance(value, str):
        return QQ(value)
    try:
        parent = value.parent()
        if parent.is_exact():
            return QQ(value)
    except (AttributeError, TypeError, ValueError):
        pass
    if isinstance(value, int):
        return QQ(value)
    raise TypeError(
        "Quaternion coordinates and Hilbert-symbol parameters must be exact; "
        "use QQ values or strings instead of floating-point numbers."
    )


def magma_iota_matrix_from_cm_data(coords, RF, cm_metadata):
    """
    Reconstruct Magma's standard Fuchsian matrix representation.

    This is for coordinates x + y*i + z*j + w*i*j in a quaternion algebra
    (a,b) with a > 0. For a nonstandard basis, pass a custom
    coordinate_matrix_function to generate_cm_image_table.
    """
    algebra_i2 = cm_metadata["algebra_i2"]
    algebra_j2 = cm_metadata["algebra_j2"]
    if algebra_i2 is None or algebra_j2 is None:
        raise ValueError(
            "cm_data must include algebra_i2 and algebra_j2, unless a custom "
            "coordinate_matrix_function is supplied."
        )

    x0, x1, x2, x3 = [_cm_exact_rational(c) for c in coords]
    qa = _cm_exact_rational(algebra_i2)
    qb = _cm_exact_rational(algebra_j2)
    if qa <= 0:
        raise ValueError(
            "The standard matrix formula requires algebra_i2 > 0; got %s. "
            "Supply coordinate_matrix_function for another presentation." % qa
        )

    sqrt_a = RF(qa).sqrt()
    sqrt_abs_b = RF(abs(qb)).sqrt()
    sign_b = RF(1) if qb > 0 else RF(-1)
    return matrix(RF, 2, 2, [
        RF(x0) + RF(x1) * sqrt_a,
        (RF(x2) + RF(x3) * sqrt_a) * sqrt_abs_b,
        sign_b * (RF(x2) - RF(x3) * sqrt_a) * sqrt_abs_b,
        RF(x0) - RF(x1) * sqrt_a,
    ])


def _cm_fixed_point_upper_half_plane(M, CC):
    """Return the fixed point of M lying in the upper half-plane."""
    aa, bb, cc, dd = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    discriminant = (dd - aa) ** 2 + 4 * bb * cc
    imaginary_unit = CC(-1).sqrt()

    if cc == 0:
        raise ValueError("The supplied elliptic matrix has c=0.")
    root1 = (
        (aa - dd) + imaginary_unit * CC(-discriminant).sqrt()
    ) / (2 * cc)
    root2 = (
        (aa - dd) - imaginary_unit * CC(-discriminant).sqrt()
    ) / (2 * cc)
    roots = [z for z in (root1, root2) if z.imag() > 0]
    if not roots:
        raise ValueError("The supplied matrix has no fixed point in H.")
    return roots[0]


def _cm_sanitize_name(value):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))
    return text.strip("_.") or "curve"


def _cm_curve_latex(E_label):
    match = re.match(r"^([0-9]+)([A-Za-z]+)([0-9]+)$", str(E_label))
    if match:
        conductor, isogeny_class, number = match.groups()
        return r"%s%s_{%s}" % (
            conductor, isogeny_class.upper(), number
        )
    return r"\mathrm{%s}" % str(E_label).replace("_", r"\_")


def _cm_default_operator_labels(keys, metadata):
    labels = {}
    ramified = metadata.get("ramified_primes")
    ramified = [] if ramified is None else list(ramified)
    for key in keys:
        text = str(key)
        suffix = text[1:] if text.startswith("w") else text
        if suffix.isdigit():
            labels[key] = suffix
        elif text == "wp" and len(ramified) >= 1:
            labels[key] = str(ramified[0])
        elif text == "wq" and len(ramified) >= 2:
            labels[key] = str(ramified[1])
        elif text == "wpq" and len(ramified) >= 2:
            labels[key] = str(Integer(ramified[0]) * Integer(ramified[1]))
        else:
            labels[key] = text
    return labels


def _cm_is_infinity(P):
    try:
        return bool(P.is_zero())
    except (AttributeError, TypeError):
        pass
    try:
        return bool(P[2] == 0)
    except (IndexError, TypeError):
        return False


def _cm_numerical_identity_from_coordinates(x, y, min_bits=80):
    r"""
    Recognize a numerical approximation to the identity of an elliptic curve.

    ``elliptic_exponential`` can return enormous affine coordinates instead of
    the exact projective point at infinity when its argument is merely
    *numerically* equal to a period.  At the identity, the standard local
    parameter ``t = -x/y`` satisfies

        x = t^-2 + O(t^-1),    y = -t^-3 + O(t^-2).

    We therefore require all three compatible estimates
    ``log2(|x|)/2``, ``log2(|y|)/3``, and ``-log2(|x/y|)`` to exceed
    ``min_bits``.  The returned bit count is a numerical diagnostic, not a
    proof that an arbitrary finite point with enormous height is the identity.
    """
    if x is None or y is None:
        return False, None

    try:
        precision = max(
            80,
            int(x.parent().precision()),
            int(y.parent().precision()),
        )
    except (AttributeError, TypeError, ValueError):
        precision = 80

    try:
        RRp = RealField(precision)
        abs_x = RRp(abs(x))
        abs_y = RRp(abs(y))
        if abs_x <= 1 or abs_y <= 1:
            return False, None

        abs_t = RRp(abs(x / y))
        if abs_t == 0:
            return True, Integer(precision)
        if abs_t >= 1:
            return False, None

        identity_bits = Integer(min(
            log(abs_x, 2) / 2,
            log(abs_y, 2) / 3,
            -log(abs_t, 2),
        ).floor())
        return identity_bits >= Integer(min_bits), identity_bits
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
        return False, None


def _cm_tail_usable_bits(stats, Nper, prec_bits):
    arithmetic_bits = max(Integer(0), Integer(prec_bits) - 50)
    if not isinstance(stats, dict):
        return arithmetic_bits
    radius = stats.get("max_u", stats.get("max_radius_seen"))
    if radius is None:
        return arithmetic_bits
    RRp = RealField(max(80, int(prec_bits)))
    radius = RRp(radius)
    if radius <= 0:
        return arithmetic_bits
    if radius >= 1:
        return Integer(0)
    tail_bits = Integer(((-Integer(Nper) * log(radius, 2)).floor()))
    tail_bits = max(Integer(0), tail_bits - 70)
    return min(arithmetic_bits, tail_bits)


def _cm_primitive_positive(R, polynomial):
    f = R(polynomial)
    if f == 0:
        return f
    content = ZZ(f.content())
    if content:
        f //= content
    if f.leading_coefficient() < 0:
        f = -f
    return f


def _cm_relative_residual_bits(f, value):
    RRp = RealField(value.parent().precision())
    absolute_value = RRp(abs(value))
    scale = RRp(0)
    for exponent, coefficient in enumerate(f.list()):
        scale += (
            abs(ZZ(coefficient))
            * max(RRp(1), absolute_value) ** exponent
        )
    residual = RRp(abs(f(value)))
    if residual == 0:
        return Integer(value.parent().precision())
    if scale == 0:
        return Integer(0)
    ratio = residual / scale
    if ratio >= 1:
        return Integer(0)
    return Integer((-log(ratio, 2)).floor())


def _cm_cluster(values, known_bits):
    if not values:
        return []
    CCp = values[0].parent()
    RRp = RealField(CCp.precision())
    comparison_bits = max(20, min(Integer(known_bits) // 3, 160))
    relative_tolerance = RRp(2) ** (-comparison_bits)
    representatives = []
    for value in values:
        if all(
            abs(value - other)
            > relative_tolerance * max(RRp(1), abs(value), abs(other))
            for other in representatives
        ):
            representatives.append(value)
    return representatives


def _cm_algdep_precisions(known_bits, cap, ladder):
    limit = min(Integer(known_bits), Integer(cap))
    values = [Integer(bits) for bits in ladder if Integer(bits) <= limit]
    if limit >= 60 and limit not in values:
        values.append(limit)
    return sorted(set(values), reverse=True)


def _cm_coefficient_height(f):
    return max([abs(ZZ(coefficient)) for coefficient in f.list()] or [ZZ(0)])


def _cm_stability_windows(known_bits, cap, ladder, radius):
    """
    Return consecutive precision windows, from highest to lowest.

    With the default radius 1, a center p is tested at p-1, p, and p+1.
    The top window is shifted down when necessary so that it stays below the
    available/capped precision.
    """
    radius = Integer(radius)
    if radius < 1:
        raise ValueError("algdep_stability_radius must be at least 1.")
    limit = min(Integer(known_bits), Integer(cap))
    centers = {
        Integer(bits) for bits in ladder
        if Integer(bits) - radius >= 60
        and Integer(bits) + radius <= limit
    }
    if limit - radius >= 60:
        centers.add(limit - radius)
    return [
        tuple(range(
            int(center - radius),
            int(center + radius) + 1,
        ))
        for center in sorted(centers, reverse=True)
    ]


def _cm_algdep_candidates_at_precision(
    selected_x,
    use_bits,
    R,
    max_degree,
    algdep_height_bound,
    algdep_min_residual_bits,
    algdep_residual_safety_bits,
):
    """Return every irreducible candidate surviving the full-value residual."""
    working_value = ComplexField(use_bits)(selected_x)
    candidates = {}
    errors = []
    for degree in range(1, int(max_degree) + 1):
        try:
            candidate = algdep(
                working_value,
                degree,
                use_bits=use_bits,
                height_bound=algdep_height_bound,
            )
        except (ValueError, ArithmeticError) as error:
            errors.append(
                "%s bits, degree %s: %s" % (use_bits, degree, error)
            )
            continue
        if candidate is None:
            continue

        candidate = _cm_primitive_positive(R, candidate)
        if candidate == 0 or not candidate.is_irreducible():
            continue
        residual_bits = _cm_relative_residual_bits(candidate, selected_x)
        required_bits = max(
            Integer(algdep_min_residual_bits),
            Integer(use_bits) - Integer(algdep_residual_safety_bits),
        )
        if residual_bits < required_bits:
            continue

        key = tuple(ZZ(coefficient) for coefficient in candidate.list())
        old = candidates.get(key)
        if old is None:
            candidates[key] = {
                "polynomial": candidate,
                "residual_bits": residual_bits,
                "required_bits": required_bits,
                "request_degrees": {Integer(degree)},
            }
        else:
            old["request_degrees"].add(Integer(degree))
            old["residual_bits"] = max(old["residual_bits"], residual_bits)
    return candidates, errors


def _cm_default_hilbert_class_field(D, h):
    """
    Construct H_D as an absolute number field for a fundamental D.

    For nonfundamental discriminants a ring-class-field constructor must be
    supplied with ``class_field_function``.
    """
    D = Integer(D)
    h = Integer(h)
    quadratic_field = QuadraticField(D, "sqrtD")
    if Integer(quadratic_field.discriminant()) != D:
        raise ValueError(
            "D=%s is not a fundamental discriminant; supply "
            "class_field_function(D,h) for its ring class field." % D
        )
    computed_h = Integer(quadratic_field.class_number())
    if computed_h != h:
        raise ValueError(
            "Magma supplied h=%s for D=%s, but Sage computes h=%s."
            % (h, D, computed_h)
        )
    relative_field = quadratic_field.hilbert_class_field("eta")
    absolute_field = relative_field.absolute_field("theta")
    expected_degree = Integer(2) * h
    if Integer(absolute_field.absolute_degree()) != expected_degree:
        raise ArithmeticError(
            "The constructed H_%s has absolute degree %s, expected %s."
            % (D, absolute_field.absolute_degree(), expected_degree)
        )
    return absolute_field


def _cm_candidate_in_class_field(polynomial, class_field):
    """Test exactly whether Q[x]/(polynomial) embeds in class_field."""
    if polynomial.degree() == 1:
        return True, Integer(1)
    candidate_field = NumberField(
        polynomial.change_ring(QQ), "candidate_root"
    )
    embeddings = candidate_field.embeddings(class_field)
    return bool(embeddings), Integer(len(embeddings))


def _cm_recognize_x(
    selected_item,
    orbit_items,
    D,
    h,
    R,
    algdep_use_bits_cap,
    algdep_use_bits_ladder,
    algdep_height_bound,
    algdep_min_residual_bits,
    algdep_stability_radius,
    algdep_residual_safety_bits,
    max_degree,
    class_field,
    class_field_error,
    require_class_field_containment,
    require_full_orbit,
    algdep_holdout_bits=30,
    algdep_min_stable_windows=2,
    algdep_window_gap_bits=20,
):
    # Consecutive-bit stability alone accepts many continued-fraction
    # approximants to transcendental numbers.  Reserve independent digits
    # of the original input, and require stability at separated windows.
    algdep_holdout_bits = Integer(algdep_holdout_bits)
    algdep_min_stable_windows = Integer(algdep_min_stable_windows)
    algdep_window_gap_bits = Integer(algdep_window_gap_bits)
    if algdep_holdout_bits < 20:
        raise ValueError("algdep_holdout_bits must be at least 20.")
    if algdep_min_stable_windows < 2:
        raise ValueError("algdep_min_stable_windows must be at least 2.")
    if algdep_window_gap_bits < 20:
        raise ValueError("algdep_window_gap_bits must be at least 20.")
    identity_count = sum(
        bool(item.get("is_infinity", False)) for item in orbit_items
    )
    if identity_count:
        if identity_count != len(orbit_items):
            return {
                "status": "inconsistent_identity_orbit",
                "polynomial": None,
                "known_bits": min(
                    Integer(item["usable_bits"]) for item in orbit_items
                ),
                "message": (
                    "%s of %s supplied CM representatives look like O. "
                    "A Galois orbit cannot mix O with finite points."
                    % (identity_count, len(orbit_items))
                ),
            }
        return {
            "status": "identity",
            "polynomial": None,
            "known_bits": min(
                Integer(item["usable_bits"]) for item in orbit_items
            ),
            "identity_bits": selected_item.get("identity_bits"),
            "identity_method": selected_item.get(
                "identity_method", "exact_or_cached"
            ),
            "message": (
                "The image is the elliptic-curve identity O; "
                "there is no finite x-coordinate to recognize."
            ),
        }

    selected_x = selected_item["x"]
    known_bits = min(
        Integer(item["usable_bits"]) for item in orbit_items
    )
    stability_windows = _cm_stability_windows(
        known_bits - algdep_holdout_bits,
        algdep_use_bits_cap,
        algdep_use_bits_ladder,
        algdep_stability_radius,
    )
    if not stability_windows:
        return {
            "status": "insufficient_tail_precision",
            "polynomial": None,
            "known_bits": known_bits,
            "message": (
                "There are not enough usable bits for a consecutive "
                "algdep stability window."
            ),
        }

    finite_x_values = [item["x"] for item in orbit_items]
    full_numeric_orbit = (
        finite_x_values + [value.conjugate() for value in finite_x_values]
    )

    candidates_by_precision = {}
    errors = []
    stable_candidates = {}
    rejected_field = []
    rejected_orbit = []
    rejected_holdout = []
    containment_cache = {}

    for window in stability_windows:
        for use_bits in window:
            if use_bits in candidates_by_precision:
                continue
            candidates, new_errors = _cm_algdep_candidates_at_precision(
                selected_x,
                Integer(use_bits),
                R,
                max_degree,
                algdep_height_bound,
                algdep_min_residual_bits,
                algdep_residual_safety_bits,
            )
            candidates_by_precision[use_bits] = candidates
            errors.extend(new_errors)

        common_keys = set(candidates_by_precision[window[0]])
        for use_bits in window[1:]:
            common_keys.intersection_update(
                candidates_by_precision[use_bits]
            )

        accepted_in_window = []
        for key in common_keys:
            common_request_degrees = set(
                candidates_by_precision[window[0]][key]["request_degrees"]
            )
            for use_bits in window[1:]:
                common_request_degrees.intersection_update(
                    candidates_by_precision[use_bits][key]["request_degrees"]
                )
            if not common_request_degrees:
                continue
            polynomial = candidates_by_precision[window[0]][key][
                "polynomial"
            ]
            stable_entry = stable_candidates.setdefault(key, {
                "polynomial": polynomial,
                "windows": [],
                "accepted_windows": [],
            })
            stable_entry["windows"].append(window)

            holdout_required_bits = max(
                Integer(algdep_min_residual_bits),
                Integer(max(window)) + algdep_holdout_bits,
            )
            selected_residual_bits = min(
                candidates_by_precision[use_bits][key]["residual_bits"]
                for use_bits in window
            )
            if selected_residual_bits < holdout_required_bits:
                rejected_holdout.append((
                    polynomial, selected_residual_bits, holdout_required_bits
                ))
                continue

            if require_class_field_containment:
                if class_field_error is not None or class_field is None:
                    rejected_field.append((
                        polynomial,
                        "class field unavailable: %s" % class_field_error,
                    ))
                    continue
                if key not in containment_cache:
                    try:
                        containment_cache[key] = (
                            _cm_candidate_in_class_field(
                                polynomial, class_field
                            )
                        )
                    except Exception as error:
                        containment_cache[key] = (False, Integer(0))
                        rejected_field.append((
                            polynomial,
                            "embedding computation failed: %s" % error,
                        ))
                        continue
                contained, embedding_count = containment_cache[key]
                if not contained:
                    rejected_field.append((
                        polynomial,
                        "Q[x]/(f) does not embed in H_%s" % D,
                    ))
                    continue
            else:
                embedding_count = None

            orbit_required_bits = holdout_required_bits
            orbit_residuals = [
                _cm_relative_residual_bits(polynomial, value)
                for value in full_numeric_orbit
            ]
            orbit_matches = sum(
                1 for bits in orbit_residuals
                if bits >= orbit_required_bits
            )
            if require_full_orbit and orbit_matches != len(
                full_numeric_orbit
            ):
                rejected_orbit.append((
                    polynomial,
                    orbit_matches,
                    len(full_numeric_orbit),
                    orbit_required_bits,
                ))
                continue

            stable_entry["accepted_windows"].append({
                "window": window,
                "request_degrees": common_request_degrees,
                "orbit_residual_bits": min(orbit_residuals),
            })
            # Keep a fixed algdep degree bound throughout all accepted
            # precisions, and require genuinely separated windows.
            chosen_windows = []
            chosen_degree = None
            for request_degree in sorted(common_request_degrees):
                eligible = [
                    entry for entry in stable_entry["accepted_windows"]
                    if request_degree in entry["request_degrees"]
                ]
                separated = []
                for entry in eligible:
                    center = (entry["window"][0] + entry["window"][-1]) / 2
                    if not separated or all(
                        abs(center - (other["window"][0] + other["window"][-1]) / 2)
                        >= algdep_window_gap_bits
                        for other in separated
                    ):
                        separated.append(entry)
                if len(separated) >= algdep_min_stable_windows:
                    chosen_windows = separated
                    chosen_degree = request_degree
                    break
            if not chosen_windows:
                continue
            stable_precisions = tuple(sorted({
                bits for entry in chosen_windows for bits in entry["window"]
            }))
            accepted_in_window.append({
                "polynomial": polynomial,
                "stable_precisions": stable_precisions,
                "stability_windows": [entry["window"] for entry in chosen_windows],
                "request_degree": chosen_degree,
                "selected_residual_bits": selected_residual_bits,
                "orbit_residual_bits": min(entry["orbit_residual_bits"] for entry in chosen_windows),
                "orbit_matches": orbit_matches,
                "embedding_count": embedding_count,
            })

        if accepted_in_window:
            best = sorted(
                accepted_in_window,
                key=lambda entry: (
                    entry["polynomial"].degree(),
                    _cm_coefficient_height(entry["polynomial"]),
                    -entry["selected_residual_bits"],
                ),
            )[0]
            return {
                "status": "recognized",
                "polynomial": best["polynomial"],
                "known_bits": known_bits,
                "algdep_use_bits": max(best["stable_precisions"]),
                "algdep_stable_precisions": best["stable_precisions"],
                "algdep_stability_windows": best["stability_windows"],
                "algdep_request_degree": best["request_degree"],
                "algdep_holdout_bits": algdep_holdout_bits,
                "algdep_support": len(best["stable_precisions"]),
                "selected_x_residual_bits":
                    best["selected_residual_bits"],
                "orbit_min_residual_bits":
                    best["orbit_residual_bits"],
                "orbit_matches": best["orbit_matches"],
                "orbit_size_with_conjugates":
                    len(full_numeric_orbit),
                "observed_distinct_x": len(
                    _cm_cluster(full_numeric_orbit, known_bits)
                ),
                "class_field_verified":
                    bool(require_class_field_containment),
                "class_field_degree": (
                    ""
                    if class_field is None
                    else Integer(class_field.absolute_degree())
                ),
                "class_field_embedding_count":
                    best["embedding_count"],
                "numerical_warning":
                    selected_item.get("warning", ""),
            }

    attempted_precisions = sorted(candidates_by_precision)
    stable_text = [
        "%s @ %s" % (
            entry["polynomial"],
            entry["windows"],
        )
        for entry in stable_candidates.values()
    ]
    if not stable_candidates:
        status = (
            "unstable_algdep"
            if any(candidates_by_precision.values())
            else "algdep_failed"
        )
        message = (
            "No polynomial was reproduced on any consecutive precision "
            "window. Different algdep answers are rejected."
        )
    elif class_field_error is not None and require_class_field_containment:
        status = "class_field_check_failed"
        message = "Could not construct H_%s: %s" % (D, class_field_error)
    elif rejected_field and not rejected_orbit and not rejected_holdout:
        status = "class_field_mismatch"
        message = (
            "Stable algdep candidates do not define subfields of H_%s. "
            "Rejected: %s" % (
                D,
                "; ".join(
                    "%s (%s)" % (polynomial, reason)
                    for polynomial, reason in rejected_field
                ),
            )
        )
    elif rejected_orbit:
        status = "orbit_mismatch"
        message = (
            "Stable candidate(s) failed the full numerical Galois-orbit "
            "test. Rejected: %s"
            % "; ".join(
                "%s (%s/%s roots at >=%s bits)"
                % values for values in rejected_orbit
            )
        )
    elif rejected_holdout:
        status = "holdout_failed"
        message = (
            "Consecutive-bit candidates failed the %s-bit holdout test "
            "on the original numerical values. No candidate also passed "
            "%s separated stability windows."
            % (algdep_holdout_bits, algdep_min_stable_windows)
        )
    else:
        status = "unstable_algdep"
        message = (
            "No candidate passed %s consecutive-precision windows with "
            "centers separated by at least %s bits at one fixed algdep "
            "degree bound."
            % (algdep_min_stable_windows, algdep_window_gap_bits)
        )

    if errors:
        message += " Last algdep error: %s" % errors[-1]
    return {
        "status": status,
        "polynomial": None,
        "known_bits": known_bits,
        "attempted_use_bits": attempted_precisions,
        "stable_candidates": stable_text,
        "class_field_verified": False,
        "message": message,
    }


def _cm_two_torsion_annotation(E, R, f):
    if f is None:
        return ""
    x = R.gen()
    division_polynomial = (
        4 * x ** 3
        + ZZ(E.b2()) * x ** 2
        + 2 * ZZ(E.b4()) * x
        + ZZ(E.b6())
    )
    return (
        "2-torsion"
        if division_polynomial.change_ring(QQ) % f.change_ring(QQ) == 0
        else ""
    )


def _cm_display_latex(row):
    if row["status"] in ("identity", "infinity"):
        return r"\mathcal{O}\;(\text{identity})"
    if row["status"] != "recognized":
        return r"\texttt{%s}" % row["status"].replace("_", r"\_")
    polynomial = row["polynomial"]
    if polynomial.degree() == 1:
        value = -QQ(polynomial[0]) / QQ(polynomial[1])
        text = "x = %s" % latex(value)
    else:
        text = latex(polynomial)
    if row["annotation"]:
        text += r"\quad(\text{%s})" % row["annotation"]
    return text


def generate_cm_image_table(
    b,
    tau0,
    S_matrices,
    AL,
    mu,
    E_label,
    Nper,
    prec_bits,
    reduction_depth,
    cm_data,
    *,
    heegner_function=None,
    atkin_lehner_keys=None,
    operator_labels=None,
    output_dir=None,
    checkpoint_file=None,
    resume=True,
    cache_tag=None,
    coordinate_matrix_function=None,
    degree_bound_function=None,
    algdep_use_bits_cap=600,
    algdep_use_bits_ladder=(80, 120, 160, 220, 300, 400, 500, 600),
    algdep_height_bound=10**40,
    algdep_min_residual_bits=40,
    algdep_stability_radius=1,
    algdep_residual_safety_bits=20,
    require_full_orbit=True,
    require_class_field_containment=True,
    class_field_function=None,
    numerical_identity_min_bits=80,
    heegner_kwargs=None,
    verbose=False,
):
    """
    Compute CM divisor images and generate algebraic-recognition tables.

    Required inputs
    ---------------
    b
        Power-series coefficient sequence.
    tau0
        Base CM/expansion point used by heegner_image_from_series.
    S_matrices
        Fuchsian side-pairing/generator matrices.
    AL
        Dictionary of Atkin-Lehner matrices. Its keys are used by default.
    mu
        Base-point parameter expected by heegner_image_from_series.
    E_label
        Sage/Cremona elliptic-curve label, for example "15a1".
    Nper, prec_bits, reduction_depth
        Numerical parameters passed unchanged to heegner_image_from_series.
    cm_data
        Dictionary made by make_magma_cm_data or
        magma_cm_data_from_namespace.

    Returns
    -------
    A dictionary containing rows, output paths, and checkpoint path.
    """
    metadata = _normalize_cm_data(cm_data)
    Nper = Integer(Nper)
    prec_bits = Integer(prec_bits)
    reduction_depth = Integer(reduction_depth)
    numerical_identity_min_bits = Integer(numerical_identity_min_bits)
    algdep_stability_radius = Integer(algdep_stability_radius)
    algdep_residual_safety_bits = Integer(algdep_residual_safety_bits)
    if numerical_identity_min_bits < 20:
        raise ValueError(
            "numerical_identity_min_bits must be at least 20."
        )
    if algdep_stability_radius < 1:
        raise ValueError("algdep_stability_radius must be at least 1.")
    if algdep_residual_safety_bits < 5:
        raise ValueError(
            "algdep_residual_safety_bits must be at least 5."
        )

    if len(b) < Nper:
        raise ValueError(
            "The coefficient sequence has length %s, but Nper=%s."
            % (len(b), Nper)
        )
    if not isinstance(AL, dict):
        raise TypeError("AL must be a dictionary of Atkin-Lehner matrices.")

    keys = (
        tuple(AL.keys())
        if atkin_lehner_keys is None
        else tuple(atkin_lehner_keys)
    )
    missing_keys = [key for key in keys if key not in AL]
    if missing_keys:
        raise ValueError(
            "AL is missing keys %s; available keys are %s."
            % (
                ", ".join(map(str, missing_keys)),
                ", ".join(map(str, AL.keys())),
            )
        )

    if heegner_function is None:
        heegner_function = globals().get("heegner_image_from_series")
    if heegner_function is None:
        raise ValueError(
            "Load heegner_image_from_series.sage first, or pass "
            "heegner_function=heegner_image_from_series."
        )

    if coordinate_matrix_function is None:
        coordinate_matrix_function = magma_iota_matrix_from_cm_data

    labels = _cm_default_operator_labels(keys, metadata)
    if operator_labels is not None:
        labels.update(operator_labels)

    curve_name = str(E_label)
    if output_dir is None:
        output_dir = "cm_results_%s" % _cm_sanitize_name(curve_name)
    os.makedirs(output_dir, exist_ok=True)
    if checkpoint_file is None:
        checkpoint_file = os.path.join(
            output_dir, "cm_image_checkpoint.sobj"
        )

    if cache_tag is None:
        cache_tag = "%s|%s|%s|%s" % (
            curve_name,
            metadata.get("algebra_i2"),
            metadata.get("algebra_j2"),
            ",".join(map(str, keys)),
        )

    raw = {}
    if resume and os.path.exists(checkpoint_file):
        loaded = load(checkpoint_file)
        if (
            isinstance(loaded, dict)
            and loaded.get("__cm_image_table_package__")
            == CM_IMAGE_TABLE_PACKAGE_VERSION
        ):
            raw = loaded.get("raw", {})
        else:
            warnings.warn(
                "Ignoring an incompatible checkpoint at %s."
                % checkpoint_file
            )

    def save_checkpoint():
        save({
            "__cm_image_table_package__":
                CM_IMAGE_TABLE_PACKAGE_VERSION,
            "raw": raw,
        }, checkpoint_file)

    # Upgrade already-computed checkpoint entries in place.  Earlier versions
    # recorded numerical approximations to O as gigantic affine coordinates,
    # which then produced the misleading status ``algdep_failed``.
    checkpoint_upgraded = False
    for item in raw.values():
        if item.get("is_infinity", False):
            if item.get("point_kind") != "identity":
                item["point_kind"] = "identity"
                item["identity_method"] = item.get(
                    "identity_method", "exact_or_cached"
                )
                checkpoint_upgraded = True
            continue
        numerical_identity, identity_bits = (
            _cm_numerical_identity_from_coordinates(
                item.get("x"),
                item.get("y"),
                numerical_identity_min_bits,
            )
        )
        if identity_bits is not None and item.get("identity_bits") != identity_bits:
            item["identity_bits"] = identity_bits
            checkpoint_upgraded = True
        if numerical_identity:
            item["is_infinity"] = True
            item["point_kind"] = "identity"
            item["identity_method"] = "local_parameter"
            checkpoint_upgraded = True
    if checkpoint_upgraded:
        save_checkpoint()

    def cache_key(record, operator):
        return (
            "cm-image-table",
            str(cache_tag),
            prec_bits,
            Nper,
            reduction_depth,
            Integer(record["D"]),
            Integer(record["orbit_index"]),
            str(operator),
        )

    CCp = ComplexField(prec_bits)
    RFp = RealField(prec_bits)
    extra_kwargs = {} if heegner_kwargs is None else dict(heegner_kwargs)

    for D in sorted(metadata["records_by_D"], key=lambda value: abs(value)):
        records = metadata["records_by_D"][D]
        for record in records:
            for operator in keys:
                key = cache_key(record, operator)
                if key in raw:
                    continue

                cm_matrix = coordinate_matrix_function(
                    record["mu_coords"], RFp, metadata
                )
                tau_D = _cm_fixed_point_upper_half_plane(cm_matrix, CCp)

                call_kwargs = {
                    "b": b,
                    "tau0": tau0,
                    "S_matrices": S_matrices,
                    "AL": AL,
                    "mu": mu,
                    "E_label": E_label,
                    "Nper": Nper,
                    "prec_bits": prec_bits,
                    "reduction_depth": reduction_depth,
                    "verbose": verbose,
                }
                call_kwargs.update(extra_kwargs)
                answer = heegner_function(
                    operator, tau_D, **call_kwargs
                )

                point = answer["P"]
                exact_infinity = _cm_is_infinity(point)
                point_x = None if exact_infinity else CCp(point[0])
                point_y = None if exact_infinity else CCp(point[1])
                numerical_identity, identity_bits = (
                    _cm_numerical_identity_from_coordinates(
                        point_x,
                        point_y,
                        numerical_identity_min_bits,
                    )
                )
                is_infinity = exact_infinity or numerical_identity
                stats = answer.get("stats", {})
                item = {
                    "D": Integer(record["D"]),
                    "h": Integer(record["h"]),
                    "orbit_index": Integer(record["orbit_index"]),
                    "operator": str(operator),
                    "prec_bits": prec_bits,
                    "Nper": Nper,
                    "reduction_depth": reduction_depth,
                    "tauD": CCp(tau_D),
                    "zE": CCp(answer["zE"]),
                    "is_infinity": is_infinity,
                    "point_kind": "identity" if is_infinity else "finite",
                    "identity_method": (
                        "exact_projective"
                        if exact_infinity
                        else (
                            "local_parameter"
                            if numerical_identity
                            else None
                        )
                    ),
                    "identity_bits": identity_bits,
                    # Keep the numerical coordinates for diagnostics and for
                    # reclassification if the threshold is changed later.
                    "x": point_x,
                    "y": point_y,
                    "stats": stats,
                    "usable_bits": _cm_tail_usable_bits(
                        stats, Nper, prec_bits
                    ),
                }
                if stats.get("bad_segments", 0):
                    item["warning"] = (
                        "Endpoint reduction reported bad_segments=%s."
                        % stats["bad_segments"]
                    )
                raw[key] = item
                save_checkpoint()
                print(
                    "computed E=%s, D=%s, class=%s/%s, %s, usable~%s bits"
                    % (
                        curve_name,
                        record["D"],
                        record["orbit_index"],
                        record["h"],
                        operator,
                        item["usable_bits"],
                    )
                )

    R = PolynomialRing(ZZ, "x")
    E = EllipticCurve(E_label)
    rows = []
    class_field_cache = {}
    class_field_error_cache = {}
    for D in sorted(metadata["records_by_D"], key=lambda value: abs(value)):
        records = metadata["records_by_D"][D]
        h = Integer(records[0]["h"])
        if require_class_field_containment and D not in class_field_cache:
            builder = (
                _cm_default_hilbert_class_field
                if class_field_function is None
                else class_field_function
            )
            try:
                class_field_cache[D] = builder(Integer(D), h)
                class_field_error_cache[D] = None
            except Exception as error:
                class_field_cache[D] = None
                class_field_error_cache[D] = str(error)
                warnings.warn(
                    "Could not construct the class field for D=%s: %s"
                    % (D, error)
                )
        class_field = class_field_cache.get(D)
        class_field_error = class_field_error_cache.get(D)
        for operator in keys:
            orbit_items = [
                raw[cache_key(record, operator)] for record in records
            ]
            max_degree = (
                Integer(2) * h
                if degree_bound_function is None
                else Integer(degree_bound_function(D, h, operator))
            )
            recognition = _cm_recognize_x(
                orbit_items[0],
                orbit_items,
                Integer(D),
                h,
                R,
                Integer(algdep_use_bits_cap),
                tuple(algdep_use_bits_ladder),
                Integer(algdep_height_bound),
                Integer(algdep_min_residual_bits),
                algdep_stability_radius,
                algdep_residual_safety_bits,
                max_degree,
                class_field,
                class_field_error,
                bool(require_class_field_containment),
                bool(require_full_orbit),
            )
            polynomial = recognition.get("polynomial")
            rows.append({
                "E": curve_name,
                "D": Integer(D),
                "h": h,
                "operator": str(operator),
                "operator_label": str(labels.get(operator, operator)),
                "status": recognition["status"],
                "polynomial": polynomial,
                "annotation": _cm_two_torsion_annotation(E, R, polynomial),
                "recognition": recognition,
            })

    csv_path = os.path.join(output_dir, "cm_x_minpolys.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "E", "D", "h", "Atkin-Lehner", "display_label", "status",
            "minimal_polynomial", "annotation", "known_bits",
            "identity_bits", "identity_method",
            "algdep_use_bits", "algdep_stable_precisions",
            "stable_candidates",
            "selected_x_residual_bits", "orbit_min_residual_bits",
            "orbit_matches", "orbit_size_with_conjugates",
            "class_field_verified", "class_field_degree",
            "class_field_embedding_count",
            "attempted_use_bits", "message",
        ])
        for row in rows:
            recognition = row["recognition"]
            writer.writerow([
                row["E"],
                row["D"],
                row["h"],
                row["operator"],
                row["operator_label"],
                row["status"],
                (
                    "" if row["polynomial"] is None
                    else str(row["polynomial"])
                ),
                row["annotation"],
                recognition.get("known_bits", ""),
                recognition.get("identity_bits", ""),
                recognition.get("identity_method", ""),
                recognition.get("algdep_use_bits", ""),
                recognition.get("algdep_stable_precisions", ""),
                recognition.get("stable_candidates", ""),
                recognition.get("selected_x_residual_bits", ""),
                recognition.get("orbit_min_residual_bits", ""),
                recognition.get("orbit_matches", ""),
                recognition.get("orbit_size_with_conjugates", ""),
                recognition.get("class_field_verified", ""),
                recognition.get("class_field_degree", ""),
                recognition.get("class_field_embedding_count", ""),
                recognition.get("attempted_use_bits", ""),
                recognition.get("message", ""),
            ])

    latex_path = os.path.join(output_dir, "cm_x_minpolys.tex")
    with open(latex_path, "w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{c r c c l}\n")
        stream.write("\\hline\n")
        stream.write(
            r"$E$ & $D$ & $h$ & $\tau$ & "
            r"minimal polynomial of $x(P(D,\tau))$ \\" + "\n"
        )
        stream.write("\\hline\n")
        for index, row in enumerate(rows):
            first_for_D = (
                index == 0 or rows[index - 1]["D"] != row["D"]
            )
            E_tex = (
                "$%s$" % _cm_curve_latex(row["E"])
                if first_for_D else ""
            )
            D_tex = "$%s$" % row["D"] if first_for_D else ""
            h_tex = "$%s$" % row["h"] if first_for_D else ""
            stream.write(
                "%s & %s & %s & $\\tau_{%s}$ & $%s$ \\\\\n"
                % (
                    E_tex,
                    D_tex,
                    h_tex,
                    row["operator_label"],
                    _cm_display_latex(row),
                )
            )
        stream.write("\\hline\n")
        stream.write("\\end{tabular}\n")

    sage_path = os.path.join(output_dir, "cm_x_minpolys.sage")
    with open(sage_path, "w", encoding="utf-8") as stream:
        stream.write(
            "# Numerical candidates; inspect CSV diagnostics before use.\n"
        )
        stream.write("R.<x> = PolynomialRing(ZZ)\n")
        stream.write("cm_x_minpolys = {\n")
        for row in rows:
            polynomial = (
                "None" if row["polynomial"] is None
                else str(row["polynomial"])
            )
            stream.write(
                "    (%s, %r): {'h': %s, 'status': %r, "
                "'polynomial': %s},\n"
                % (
                    row["D"],
                    row["operator"],
                    row["h"],
                    row["status"],
                    polynomial,
                )
            )
        stream.write("}\n")

    rows_path = os.path.join(output_dir, "cm_x_minpolys_rows.sobj")
    save(rows, rows_path)

    recognized = sum(row["status"] == "recognized" for row in rows)
    identity = sum(
        row["status"] in ("identity", "infinity") for row in rows
    )
    unresolved = [
        row for row in rows
        if row["status"] not in ("recognized", "identity", "infinity")
    ]
    print(
        "Finished %s: %s recognized, %s identity points, %s unresolved."
        % (curve_name, recognized, identity, len(unresolved))
    )
    print("Wrote %s" % output_dir)
    for row in unresolved:
        message = row["recognition"].get("message", "")
        print(
            "  D=%s, %s: %s%s"
            % (
                row["D"],
                row["operator"],
                row["status"],
                "" if not message else " -- " + message,
            )
        )

    return {
        "rows": rows,
        "paths": {
            "csv": csv_path,
            "latex": latex_path,
            "sage": sage_path,
            "rows": rows_path,
        },
        "checkpoint": checkpoint_file,
        "metadata": metadata,
    }


def _cm_absolute_degree(field):
    if field is QQ:
        return Integer(1)
    try:
        return Integer(field.absolute_degree())
    except AttributeError:
        return Integer(field.degree())


def _cm_field_polynomial(field):
    if field is QQ:
        return None
    try:
        return field.defining_polynomial()
    except AttributeError:
        return field.polynomial()


def _cm_subfield_description(H, elements, name):
    """
    Return the exact subfield generated by ``elements`` inside H.

    Sage returns the descended elements and the inclusion map into H, so the
    equalities with the original ambient elements can be checked exactly.
    """
    field, descended, inclusion = H.subfield_from_elements(
        list(elements),
        name=name,
        polred=True,
    )
    if len(descended) != len(elements):
        raise ArithmeticError(
            "Subfield descent returned the wrong number of elements."
        )
    if any(
        inclusion(descended[index]) != H(elements[index])
        for index in range(len(elements))
    ):
        raise ArithmeticError(
            "A descended coordinate does not map back to its element in H_D."
        )
    return {
        "field": field,
        "degree": _cm_absolute_degree(field),
        "defining_polynomial": _cm_field_polynomial(field),
        "elements": list(descended),
        "inclusion": inclusion,
    }


def _cm_safe_field_name(prefix, D, operator, index):
    text = "%s_%s_%s_%s" % (
        prefix,
        "m%s" % abs(Integer(D)) if Integer(D) < 0 else Integer(D),
        _cm_sanitize_name(operator),
        Integer(index),
    )
    return re.sub(r"[^A-Za-z0-9_]", "_", str(text))


def _cm_recover_row_points(row, H):
    """
    Recover all exact points associated with one recognized x-polynomial.

    The output contains both signs of y and all Galois-conjugate x-roots.
    """
    E = EllipticCurve(str(row["E"]))
    EH = E.change_ring(H)
    polynomial = row.get("polynomial")
    if polynomial is None:
        raise ValueError("A recognized row has no x-polynomial.")

    RX = PolynomialRing(H, "X")
    x_roots = RX(polynomial).roots(multiplicities=False)
    expected_x_roots = Integer(polynomial.degree())
    if len(x_roots) != expected_x_roots:
        raise ValueError(
            "x-field error: f has %s roots in H_%s, expected %s."
            % (len(x_roots), row["D"], expected_x_roots)
        )

    RY = PolynomialRing(H, "Y")
    Y = RY.gen()
    recovered = []
    point_index = Integer(0)

    for xH in x_roots:
        rhs = (
            xH ** 3
            + H(E.a2()) * xH ** 2
            + H(E.a4()) * xH
            + H(E.a6())
        )
        linear_y = H(E.a1()) * xH + H(E.a3())
        y_polynomial = Y ** 2 + linear_y * Y - rhs
        y_roots = y_polynomial.roots(multiplicities=False)
        if not y_roots:
            raise ValueError(
                "y-field error: the Weierstrass equation has no y-root "
                "in H_%s above x=%s." % (row["D"], xH)
            )

        for yH in y_roots:
            point_index += 1
            point_H = EH([xH, yH, H(1)])

            name_x = _cm_safe_field_name(
                "x", row["D"], row["operator"], point_index
            )
            name_y = _cm_safe_field_name(
                "y", row["D"], row["operator"], point_index
            )
            name_point = _cm_safe_field_name(
                "p", row["D"], row["operator"], point_index
            )

            x_field = _cm_subfield_description(H, [xH], name_x)
            y_field = _cm_subfield_description(H, [yH], name_y)
            point_field = _cm_subfield_description(
                H, [xH, yH], name_point
            )
            L = point_field["field"]
            xL, yL = point_field["elements"]
            EL = E.change_ring(L)
            point_L = EL([xL, yL, L(1)])
            inclusion = point_field["inclusion"]

            mapped_point = EH([
                inclusion(point_L[0]),
                inclusion(point_L[1]),
                H(1),
            ])
            if mapped_point != point_H:
                raise ArithmeticError(
                    "Point-field error: the descended point does not map "
                    "back to the point over H_D."
                )
            if _cm_absolute_degree(H) % point_field["degree"] != 0:
                raise ArithmeticError(
                    "Point-field error: Q(x,y) has degree %s, which does "
                    "not divide [H_D:Q]=%s."
                    % (
                        point_field["degree"],
                        _cm_absolute_degree(H),
                    )
                )

            recovered.append({
                "point_index": point_index,
                "ambient_point": point_H,
                "point": point_L,
                "curve": EL,
                "x": xL,
                "y": yL,
                "x_field": x_field,
                "y_field": y_field,
                "point_field": point_field,
            })

    return {
        "E": str(row["E"]),
        "D": Integer(row["D"]),
        "h": Integer(row["h"]),
        "operator": str(row["operator"]),
        "status": "algebraic_points_recovered",
        "ambient_class_field": H,
        "ambient_class_field_degree": _cm_absolute_degree(H),
        "x_polynomial": polynomial,
        "points": recovered,
        "message": "",
    }


def recover_recognized_points_over_subfields(
    result,
    *,
    save_path=None,
    verbose=True,
):
    r"""
    Recover exact points over their actual coordinate fields.

    For every row of ``result`` with status ``recognized``, this function:

    1. constructs the exact Hilbert class field H_D;
    2. factors the recognized x-polynomial in H_D;
    3. solves the Weierstrass equation for y in H_D;
    4. computes Q(x), Q(y), and Q(x,y) with exact inclusion maps into H_D;
    5. constructs P on E/Q(x,y) and verifies its image in E/H_D.

    A failure is returned as a row with status beginning ``error_`` instead
    of being silently accepted.  The output may be saved as a Sage ``.sobj``
    file with ``save_path``.
    """
    if not isinstance(result, dict) or "rows" not in result:
        raise TypeError(
            "result must be the dictionary returned by "
            "generate_cm_image_table(...)."
        )

    recognized_rows = [
        row for row in result["rows"]
        if row.get("status") == "recognized"
    ]
    class_fields = {}
    output_rows = []

    for row in recognized_rows:
        D = Integer(row["D"])
        h = Integer(row["h"])
        key = (D, h)

        if key not in class_fields:
            try:
                class_fields[key] = (
                    _cm_default_hilbert_class_field(D, h), None
                )
            except Exception as error:
                class_fields[key] = (None, str(error))

        H, class_field_error = class_fields[key]
        if class_field_error is not None:
            recovered_row = {
                "E": str(row["E"]),
                "D": D,
                "h": h,
                "operator": str(row["operator"]),
                "status": "error_class_field",
                "points": [],
                "message": class_field_error,
            }
        else:
            try:
                recovered_row = _cm_recover_row_points(row, H)
            except Exception as error:
                message = str(error)
                if message.startswith("x-field error"):
                    status = "error_x_not_in_class_field"
                elif message.startswith("y-field error"):
                    status = "error_y_not_in_class_field"
                elif message.startswith("Point-field error"):
                    status = "error_point_field_descent"
                else:
                    status = "error_algebraic_recovery"
                recovered_row = {
                    "E": str(row["E"]),
                    "D": D,
                    "h": h,
                    "operator": str(row["operator"]),
                    "status": status,
                    "ambient_class_field": H,
                    "ambient_class_field_degree":
                        _cm_absolute_degree(H),
                    "x_polynomial": row.get("polynomial"),
                    "points": [],
                    "message": message,
                }

        output_rows.append(recovered_row)
        if verbose:
            print(
                "E=%s, D=%s, %s: %s"
                % (
                    recovered_row["E"],
                    recovered_row["D"],
                    recovered_row["operator"],
                    recovered_row["status"],
                )
            )
            for point_data in recovered_row.get("points", []):
                print(
                    "  point %s: deg Q(x)=%s, deg Q(y)=%s, "
                    "deg Q(x,y)=%s"
                    % (
                        point_data["point_index"],
                        point_data["x_field"]["degree"],
                        point_data["y_field"]["degree"],
                        point_data["point_field"]["degree"],
                    )
                )
                print("    Q(x)   =", point_data["x_field"]["field"])
                print("    Q(y)   =", point_data["y_field"]["field"])
                print(
                    "    Q(x,y) =",
                    point_data["point_field"]["field"],
                )
                print("    point =", point_data["point"])
            if recovered_row.get("message"):
                print("  " + recovered_row["message"])

    output = {
        "rows": output_rows,
        "recognized_row_count": Integer(len(recognized_rows)),
        "recovered_row_count": Integer(sum(
            row["status"] == "algebraic_points_recovered"
            for row in output_rows
        )),
        "error_row_count": Integer(sum(
            row["status"] != "algebraic_points_recovered"
            for row in output_rows
        )),
    }
    if save_path is not None:
        save(output, save_path)
        output["save_path"] = save_path
        if verbose:
            print("Saved exact algebraic point data to %s" % save_path)
    return output
