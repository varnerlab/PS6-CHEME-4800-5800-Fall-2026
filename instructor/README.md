# PS6 instructor notes

The student task is fixed: simulate the regulated feedback (memory) motif of Alon
by solving the full nonlinear mRNA and protein balances with an ODE solver
(`Tsit5`), using the course's gene-network model. The network and its parameters
come from a TOML file in the course format. Parts 1–2 require three functions
(`control`, `balances!`, `simulate`) and two written responses, graded out of four
points using 40 tests. There is no Advanced track. The loader (`load_network`),
`bound_fraction`, `inducer`, and the report code (`src/Reporting.jl`) are
supplied. The design is recorded in
[plans/2026-11-21-ps6-ode-redesign.md](plans/2026-11-21-ps6-ode-redesign.md);
[DESIGN.md](DESIGN.md) is the earlier design, with linear translation and the
bilinear recursion, and is stale. Both are local files, ignored by Git, because
they list the ideas expected in the answers.

## Reviewing the assignment

- Read [the problem statement](../README.md), [rubric](../RUBRIC.md), and
  [questions](../responses.md).
- Inspect the reference reports in [reference/](reference/):
  [memory.png](reference/memory.png) and the CSV tables.
- The local, Git-ignored [reference solution](../solution/README.md) contains
  completed code, sample written responses, and the parameter tuning record.
- [The student ZIP](artifacts/PS6-CHEME-4800-5800-Fall-2026-student.zip)
  contains starter code and data, with no instructor files, outputs, or solutions.

Release is Saturday, November 21, 2026. The initial submission deadline is
11:59 PM ET on Saturday, December 5, 2026. No commit, push, tag, GitHub release,
or Canvas upload is performed by these scripts.

Students download the **Source code (zip)** of the tagged GitHub release
`ps6-cheme-4800-5800-2026.1` of `varnerlab/PS6-CHEME-4800-5800-Fall-2026`;
`.gitattributes` keeps instructor material out of that archive. Pushing the tag
runs the [release Action](../.github/workflows/release-ps6.yml), which runs
[validate_release.py](validate_release.py) on the source archive and creates a
**draft** release with the [release notes](release-notes.md). That validator
needs no reference solution; run it locally with
`python3 instructor/validate_release.py` before tagging (it needs a committed
Git tree). Review and publish the draft, then post the
[Canvas description](canvas-assignment-description.html), which links to the
same tag. Later corrections use a new tag.

## Reproducing the checks

To run the reference solution directly from the assignment root:

```bash
julia --startup-file=no --project=. check_submission.jl --solution
```

This runs all 40 tests, checks `solution/responses.md`, and writes the reports to
`solution/outputs/` and a manifest to `solution/MANIFEST.txt`. The student files
and outputs stay unchanged. The same flag works with `testme_part_1.jl`,
`testme_part_2.jl`, and `runmotif.jl`. The default remains the student
implementation; a missing local solution causes an error rather than falling back
to student code.

From the assignment root, install the Julia dependencies and run:

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
python3 instructor/validate.py
python3 instructor/build_student_package.py
```

The validator first runs [check_versions.jl](check_versions.jl), which checks
that every package shared with the course environment has the course version; it
skips if the course repository is not next to this one. It also runs the local,
Git-ignored checks that the solution README describes, when they are present. The
validator then runs [test_supplied.jl](test_supplied.jl), separate from student
scores, once without and once with `--solution`. It checks the loader (the fields
of `memory.toml`, the kinetic limits and time constants, the weights from `dG`,
the loss constants, unknown regulators, missing fields, values that must be
positive, the temperature, other inducer names, the `[inputs.I]` table, Windows
line endings, and the no-loop file), the helpers (bound fractions, negative
concentrations, and the inducer pulses), and, with `--solution`, the three steady
states. The validator then copies the student files into temporary directories
and substitutes the reference implementation only in those copies. It checks all
40 tests, the driver, and the
checker's handling of unfinished code, partial work, syntax errors, missing
responses, Windows line endings, runs from another folder, and reporting errors.
Logs are in `instructor/validation-output/`; reference reports are copied to
`instructor/reference/`. Both directories are ignored by Git and excluded from
the student package.

## Reference results

The reference numbers are in the Git-ignored
[reference answers](../solution/responses.md) and
[solution README](../solution/README.md). They are the numbers students report,
so they stay off GitHub until the reference solution is released.

## Modeling choices

- **Units.** Intracellular concentrations, `G`, `K_X`, `K_L`, and `K` are in nmol/gDW;
  `Vmax_X` and `Vmax_L` are in nmol/gDW/min; the inducer and `K_I` are in µM, in
  the medium; `dG` is in kJ/mol; the temperature is in K; time is in minutes.
- **Model (instructor, November 21, 2026).** The course's gene-network balances:
  transcription is `Vmax_X·G/(τ_X·K_X + (1 + τ_X)·G + O_X)` times the control variable,
  translation is `Vmax_L·m/(τ_L·K_L + (1 + τ_L)·m + O_L)`, where `O_X` and `O_L` are the
  competition of the other genes for RNA polymerase and of the other mRNAs for ribosomes,
  and mRNA and protein are lost by degradation and dilution. The loader computes `G`,
  `Vmax`, `τ`, and `K` from the cell, RNA polymerase, and ribosome values in `[parameters]`.
  There is no translational control. The control variable follows the active-state rule:
  the empty promoter has weight 1, and each listed state has weight
  `W = exp(−dG/RT)` times the bound fraction of each regulator bound in it.
- **Values.** The cell, RNA polymerase, ribosome, and gene values, the mRNA half-life,
  the doubling time, the temperature, the leak and activated-state `dG`, and `K_I` are
  the course values. The protein half-life is 35 min, for a growing cell; the course's
  70 h is a cell-free value. The X-bound state's
  `dG`, K(X), n(X), and the partner K were tuned with the Git-ignored
  `instructor/tune_parameters.jl`; the tuning record is in the solution README,
  because it gives the answers.
- **Solver.** `simulate` uses `Tsit5` with `reltol = 1e-8`, `abstol = 1e-12`,
  `saveat = 1.0`, and `tstops` at the pulse edges inside the run. The tests compare
  with stored values at `rtol = 1e-6`. Checked by hand: a `simulate` without
  `tstops` fails the half-minute pulse case, default tolerances fail 9 cases, and a
  control variable read only at `x0` fails 3 cases.
- **Remembered.** A pulse is remembered if protein Y is above the middle steady
  state at 900 min. The middle steady state is found on the line `m_Y = m_Z`,
  `p_Y = p_Z` by a scan over `m_Y` and bisection in the supplied code, which gets every
  rate from the student's `balances!`.
- **Exact numbers.** The tests compare with literals stored in `test_support.jl`,
  so no reference code ships. The Git-ignored `instructor/make_stored_values.jl`
  prints those literals from the reference solution.
- **Verified facts.** Alon (2019), *An Introduction to Systems Biology*, 2nd
  edition, chapter 5, checked against the Taylor & Francis listing. Alon (2007),
  *Nature Reviews Genetics* 8(6), 450–461, DOI 10.1038/nrg2102: volume, issue, and
  DOI checked; the motif text was not fetched (the server refused).

## Preparing a release

Review the dates in the main README and the student workload, then run the
validation and build scripts. The builder refuses a release package with dates
still marked to be announced. It uses an explicit allowlist and checks that the
three starter errors remain in `src/Compute.jl`.

The reference code, sample answers, and generated reference outputs are kept
outside the student archive.
