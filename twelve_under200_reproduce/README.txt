SOURCE-ONLY REPRODUCTION: TWELVE CM-DIVISOR TABLES, |D|<200
=======================================================

This package contains only the code, exact CM/geometry inputs, and table
renderer needed for the computation. No old numerical/result SOBJ files,
tables, logs, historical audits, notebooks or baseline bundles are included.
The generated checkpoints are OUTPUTS of your run, not required inputs.

Scope
-----
Targets:14a1,15a1,21a1,26a1,26b1,33a1,34a1,38a1,38b1,46a1,58a1,58b1.
Negative fundamental discriminants |D|<200; no class-number cutoff; all
ideal classes and local orientations; all three nonidentity AL operators.
Each image is that of [w*tau_D]-[tau_D] on the Sage X_0(N)-optimal target.
The reference table has5616 images (4268 finite,1348 identities),720 rows.

Run in Sage/Ubuntu
------------------
Extract the directory intact. In a terminal:

    cd /path/to/twelve_under200_reproduce
    sage

Then in Sage (or a Sage Jupyter notebook whose cwd is this directory):

    load("reproduce.sage")
    results = run_all()

This is a FULL from-source computation, not the earlier incremental job.
Coefficient generation and algebraic recovery can be slow. Files are
checkpointed under results/; rerunning the same command resumes completed
work with matching input/settings fingerprints. Loading a file by itself
does not start the full computation. Do not mix these generated caches
with files from the previous large bundle.

For a single target first:

    results = run_all(labels=("15a1",))

To choose another output directory or override a pair:

    results = run_all(labels=("26b1",), output_dir="my_results",
                      pass_specs={"26b1": (3000,3200,1600)})

Each tuple is (lower Nper, higher Nper, working bits). Defaults are in
CM_PASS_SPECS in reproduce.sage. They use the strongest final settings
uniformly over the full range, rather than depending on old-range caches.
In particular26b1 uses1600 bits; its large-height cases cannot be addressed
by merely padding1000-bit numbers. The scale constants and coefficients
are genuinely regenerated at guarded higher precision.

Generated tables and exact points
--------------------------------
run_all() automatically writes into results/:
  cm_divisors_twelve_under200.csv             compact table
  cm_divisors_twelve_under200_all_points.csv  all individual images
  cm_divisors_twelve_under200_data.json       exact coefficients/metadata
  cm_divisors_twelve_under200.tex             editable LaTeX table
  cm_divisors_twelve_under200_validation.json validation summary
  xN_cm_images_LABEL_result.sobj             exact points/fields/maps
alongside its newly generated coefficient/period/numerical checkpoints.
A subset uses the filename prefix cm_divisors_subset instead.

To regenerate tables from this package's completed result files only:

    load("reproduce.sage")
    load("export_table.sage")
    table = export_table()

The exporter checks full class/orientation coverage, exact point/field
inclusions, both coordinate matches, independent truncation consistency,
and the exact Atkin-Lehner point identities. Its built-in SHA256 fingerprint
compares all720 primitive x-polynomials/statuses with the published table,
without loading any old output file. A mismatch is reported, not ignored.

To retrieve a point and its minimal coordinate fields after the run:

    entry = cm_divisor_point("26b1", -187, "w13",
                              orientation_index=1, class_index=1)
    print(entry["point_field"]["field"], entry["point"])
    print(entry["x_field"]["field"], entry["y_field"]["field"])
    inclusion = entry["point_field"]["inclusion"]  # Q(x,y) -> H_D
    small = compact_point(entry)                  # optional presentation change

The identity case returns status="identity", point=E(0), field=QQ.
The x/y/point fields, inclusions into H_D, numerical branch and sign of y
are retained in the Sage result, not merely printed as polynomial strings.

Runtime
-------
The original source was prepared with SageMath 10.9. This certificate package
also runs the small sanity check with SageMath 9.5, using equivalent exact
discriminant tests and quaternion inverses. The optional PDF-only renderer
is retained with the computation and exact data exporters; it is not run by
the certificate entrypoint and is not needed for the numerical checks.
See the parent README.md and results/cm_sanity.json for the current run.

Magma
-----
cm_data/xN_cm_under200.m is a standalone print-only Magma program for
each N=14,15,21,26,33,34,38,46,58. Paste one into the online calculator:
https://magma.maths.usyd.edu.au/calc/
It covers the complete range3..199, with no class-number cutoff, and
prints Sage-readable exact coordinates/order/Atkin-Lehner data. No Magma
file writing is required. Copy the printed output into xN_cm_data.sage.

The supplied cm_data/xN_cm_data.sage files already contain the exact
representatives used for the actual table, so Magma is optional for a
direct reproduction. Fresh Magma runs can choose equivalent representatives
or different class indices; use a new results directory after replacing
CM inputs. See cm_data/README.txt for the precise data contract.

Numerical method and acceptance checks
-------------------------------------
The supplied recurrences generate b_n; exact side geometry builds the
Dirichlet reducer; full side periods determine the scalar to the optimal
elliptic curve, including the native-to-optimal isogeny factor. Two different
series truncations estimate usable precision. The genus-two runner builds
one guarded high-order series and truncates its prefix for the lower pass;
the two integrations have different truncation lengths.

26b1,38b1,58b1 use checked kernel-w2 preconditioning when it improves the
reduced endpoint radii.33a1 uses w11 preconditioning for P_w3 followed by
P_w33=P_w3+(0,0). The rational2-torsion translation cancels in divisor
differences. Exact normalizer/product checks and numerical differential
eigenvalue checks are regenerated; direct-coordinate bridges are retained.
Origin and torsion cases are handled explicitly, including in the old range.

Recognition requires separated windows of consecutive algdep BIT precisions,
at least30 held-out bits, exact lifting to E(H_D), matching BOTH x and y,
and descent to Q(x),Q(y),Q(x,y). Large polynomial heights do not waive these
checks. Algebraic equations and field inclusions are exact; identification
with the analytic CM images is numerical, not interval-certified proof.

Small sanity check; not a full rerun
-----------------------------------
    load("sanity_check.sage")
    sanity_check()

This parses the source modules, checks exact CM coverage/order membership,
computes only20 coefficients per form, checks CM lifts and tiny AL/torsion
identities. It does not integrate the divisor table or run recognition.
These light checks passed during preparation. The simplified period matchers
also reproduced the saved complex normalizations in a read-only comparison.
The complete from-scratch pipeline was NOT rerun for this code-only request.

The only subdirectories are common/ (shared engines), cm_data/ (exact inputs
and optional Magma), six_curves/ and genus_two/ (recurrences/geometry/drivers).
All runtime paths are relative to this extracted package. In a notebook,
use os.chdir('/path/to/twelve_under200_reproduce') before loading.
For very large all-points CSV fields, set csv.field_size_limit(10**9) in Python.
