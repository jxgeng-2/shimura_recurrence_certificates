# Computational certificates

For **Recurrence formulas for Taylor expansions at CM points on genus-one and genus-two Shimura curves**, by Jingxuan Geng. Section, appendix and table numbers follow
`recurrence_formulas_genus_one_refined_edited.pdf` (the version checked on
2026-11-19). Appendix A contains both initial values and complete CM expansions; Appendix B contains
certification of CM coordinates.

| Requested computation | Open this file |
|---|---|
| 1. All Table 3 families | [Section_3_3 Verify Families of R(u) in Table 3.ipynb](<Section_3_3 Verify Families of R(u) in Table 3.ipynb>) |
| 2. Sections 4.1–4.2: Tu's Schwarzian formulas and recovery from coefficients and the curve | [Sections_4_1_and_4_2 Recover R(u) and Recurrence Parameters.ipynb](<Sections_4_1_and_4_2 Recover R(u) and Recurrence Parameters.ipynb>) |
| 3. Section 4.3: integrality and every row of Table 4 | [Section_4_3 Verify Integrality of c_n and Table 4.ipynb](<Section_4_3 Verify Integrality of c_n and Table 4.ipynb>) |
| 4. Section 5: the map Phi and CM divisor images | [HTML table](<Section_5 Compute Phi and CM Images/Section_5 Table of All CM Images.html>); [Sage computation](<Section_5 Compute Phi and CM Images/Compute Phi and CM Images.sage>) |
| 5. Appendix A: initial coefficients by Nelson's method | [Appendix_A Initial Values by Nelson Method/Appendix_A Compute Initial Coefficients by Nelson Method.ipynb](<Appendix_A Initial Values by Nelson Method/Appendix_A Compute Initial Coefficients by Nelson Method.ipynb>) |
| 6. Appendix A: Tau, Theta and the first 20 coefficients for every form | [Appendix_A Compute Tau Theta and First 20 Coefficients.ipynb](<Appendix_A Compute Tau Theta and First 20 Coefficients.ipynb>) |

## Run the notebooks

The notebooks were executed using SageMath 9.5 on Ubuntu 22.04 under WSL.
If Sage is already installed, keep that installation. For a new Ubuntu 22.04
installation, the relevant packages are:

```sh
sudo apt-get update
sudo apt-get install sagemath sagemath-database-elliptic-curves jupyter-notebook
```

In an **Ubuntu/WSL terminal**, enter the package directory and register the
Sage notebook kernel once. This also repairs an older kernel entry that
incorrectly launches a command named `python`:

```sh
cd /mnt/c/Users/Jingx/OneDrive/Documents/shimura_recurrence_certificates
/usr/bin/sage -python -c 'from sage.repl.ipython_kernel.install import SageKernelSpec; import site; SageKernelSpec.update(prefix=site.USER_BASE)'
/usr/bin/sage -n jupyter
```

Open a notebook, select the **SageMath** kernel, then choose **Kernel → Restart
& Run All**. In the checked run, the first three notebooks took seconds,
the Appendix A expansion notebook took about 6 minutes, and the Nelson notebook took about 12 minutes.
The explicit `/usr/bin/sage` selects Ubuntu's installed
Sage instead of a different source checkout earlier on PATH. On another
installation, substitute the path to its working Sage executable. The clean
installation recipe itself has not been tested in a newly installed system.

In the Appendix A expansion notebook, `first_coefficients("46a", N)` returns the exact normalized
coefficients `c_0,...,c_(N-1)`. They occur in
`f(z)=(1-w)^2 sum(c_n*(Theta*w)^n/n!)`. All twelve lists with `N=20` are displayed.
The notebook also checks the explicit maximal-order bases defining Gamma and
computes tau from the CM fixed-point equation. No Magma run is needed for this.

## Section 5: view or reproduce the table

Open [Section_5 Table of All CM Images.html](<Section_5 Compute Phi and CM Images/Section_5 Table of All CM Images.html>)
in a browser. This self-contained file works offline and prints directly. Its
five columns list 720 rows: rational values, minimal polynomials, points at
infinity and 2-torsion annotations. Repeated entries are listed once.

From the package directory, in an Ubuntu/WSL terminal:

```sh
cd "Section_5 Compute Phi and CM Images"
/usr/bin/sage "Compute Phi and CM Images.sage" --sanity
```

This runs small checks for all twelve forms and the three Table 5 polynomials:
`14a1, D=-163, w2`; `15a1, D=-127, w3`; and **`58a1, D=-131, w29`**.
It does not recompute their analytic CM images. The last polynomial is
`83521*x^5 - 156927*x^4 + 139246*x^3 - 185426*x^2 + 211985*x - 88687`.

To compute the entire table from the supplied inputs:

```sh
/usr/bin/sage "Compute Phi and CM Images.sage" --full
```

This computes and validates all twelve forms, saves intermediate results in
`results/`, and writes `Section_5 Table of All CM Images.html`. The full numerical
calculation is expensive; compatible checkpoints let an interrupted run resume.
No downloaded results archive or Magma installation is needed.

If `results/` already contains a completed run, rebuild the HTML without
repeating the integrations and recognition:

```sh
/usr/bin/sage "Compute Phi and CM Images.sage" --table
```

This rechecks the saved results, confirms all three current Table 5 examples,
and writes the same table. Its columns omit orientation labels; the code checks
that the repeated polynomials agree before combining them.

For a single target, start Sage in this folder and use:

```sage
load("Compute Phi and CM Images.sage")
results = run_all(labels=("58a1",))
point = cm_divisor_point("58a1", -131, "w29")
point["point"]
```

The complete run on 2026-11-19, using SageMath 9.5, covers **5,616 divisor images**
for every admissible negative fundamental discriminant with `|D| < 200`, all
twelve forms, ideal classes and three Atkin–Lehner operators. All coverage and
validation checks passed, with no unresolved entries. The minimum observed
agreement in both coordinates is 169 bits, using
`|z_num-z_exact| / max(1, |z_num|, |z_exact|)` for finite points. This is a numerical
agreement measure. The current `58a1` example was checked against that saved run.

For distribution, keep the HTML, the six `.sage` files, the fundamental-domain
JSON input, and `Generate CM Embedding Data.m`. The last file records optional
regeneration of the supplied CM inputs in Magma; it is not needed for the Sage
commands. These nine files give the table and its reproducible computation.
The large `Full CM Image Results.zip`, unpacked `results/`, logs and generated
`.sage.py` files are local working material, excluded from the sharing ZIP and
ignored by Git. Retain the large archive separately if you want the original
checkpoints and number-field maps; readers do not need it to view or recompute
the table.

## Meaning of the checks

The family and recovery notebooks identify their geometric and initial-value
assumptions. The integrality notebook verifies the polynomial-system hypotheses,
including that `Q(0)` is a unit, and gives the argument for **all** coefficients;
its conclusion is not inferred from a finite coefficient sample. Minimal prime
sets are checked using the reduced denominators of `c_2` and `c_3`.

Nelson evaluations use cutoffs 144 and 196 and report observed agreement with
the stated values. They do not provide interval bounds for the omitted
infinite lattice tail or accumulated rounding. Recurrence coefficients are
exact consequences of the stated input data. Section 5's analytic identification
remains numerical; its saved full results include two-truncation comparisons
and exact algebraic checks.

## Reuse and AI assistance

Copyright (c) 2026 Jingxuan Geng. The code and accompanying documentation are
licensed under **GPL-3.0-only**; see [LICENSE](<LICENSE>). Cited publications and
external software retain their own licenses. Please cite the paper and the
software commit or release used. SageMath supplies the exact algebra,
elliptic-curve database, power series and multiprecision arithmetic; the Nelson
notebook identifies the source of its numerical formula.

The author reports using **ChatGPT and Claude throughout work on the paper**.
Codex assisted with implementation, review, simplification, testing and
documentation of these certificates.

To share this version, commit the notebooks, the nine Section 5 files listed
above, this README, LICENSE and `.gitignore`. Saved notebook outputs belong with
the code. The sharing ZIP omits the large results archive and generated working
files.
