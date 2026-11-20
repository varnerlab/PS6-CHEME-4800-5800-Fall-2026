# PS6 instructor notes

The student task is fixed: simulate the regulated feedback (memory) motif of Alon
with the bilinear recursion. The network and its parameters come from a TOML
file. Parts 1–2 require three functions (`control`, `discretize`, `simulate`) and
two written responses, graded out of four points using 40 tests. There is no
Advanced track. The loader (`load_network`), `state_matrices`, `bound_fraction`,
`inducer`, and the report code (`src/Reporting.jl`) are supplied. The design is
recorded in [DESIGN.md](DESIGN.md) and the build plan in [plans/](plans/); both
are local files, ignored by Git, because DESIGN.md lists the ideas expected in the
answers.

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

The validator first runs [test_supplied.jl](test_supplied.jl), separate from
student scores, once without and once with `--solution`. It checks the loader
(the fields of `memory.toml`, the loss constants, unknown regulators, missing
fields, other inducer names, Windows line endings, and the no-loop file), the
helpers (`A` and `B`, bound fractions, negative levels, and the inducer pulses),
and, with `--solution`, the three steady states. The validator then
copies the student files into temporary directories and substitutes the reference
implementation only in those copies. It checks all 40 tests, the driver, and the
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

- **Units.** Levels inside the cell, `Vmax_X`, and `K` are in nmol/gDW; the
  inducer and `K_I` are in µM, in the medium; time is in minutes (instructor,
  November 19, 2026).
- **Model.** Transcription is `Vmax_X` times the control; translation is `k_L`
  times the mRNA; there is no translational control, so `A` and `B` are constant
  and all of the regulation is in the control. The control follows the
  active-state rule: the bare promoter has weight 1, and each listed state has
  weight `W` times the bound fraction of each regulator bound in it.
- **Values.** `Vmax_X`, `k_L`, the half-lives, and the doubling time are derived
  from course cell constants; [data/README.md](../data/README.md) gives the
  derivation. W₁, K(X), n(X), and K(Y, Z) were tuned with the Git-ignored
  `instructor/tune_parameters.jl`; the tuning record is in the solution README,
  because it gives the answers.
- **Time step (instructor, November 20, 2026).** Every report run uses 0.1 min,
  and the README states it. To keep PS6 short, there is no question about the
  time step and no step-size table. The tests still call `simulate` at steps of
  0.1 to 15 min, because they check the recursion, not its accuracy.
- **Remembered.** A pulse is remembered if protein Y is above the middle steady
  state at 900 min. The middle steady state is found on the line `p_Y = p_Z` by a
  scan and bisection in the supplied code.
- **Exact numbers.** The tests compare with literals stored in `test_support.jl`,
  so no reference code ships. The Git-ignored `instructor/make_stored_values.jl`
  prints those literals from the reference solution.
- **Verified facts.** Alon (2019), *An Introduction to Systems Biology*, 2nd
  edition, chapter 5, checked against the Taylor & Francis listing. Alon (2007),
  *Nature Reviews Genetics* 8(6), 450–461, DOI 10.1038/nrg2102: volume, issue, and
  DOI checked; the motif text was not fetched (the server refused). Gu, Goel, and
  Ré (2022), ICLR, arXiv 2111.00396.

## Preparing a release

Review the dates in the main README and the student workload, then run the
validation and build scripts. The builder refuses a release package with dates
still marked to be announced. It uses an explicit allowlist and checks that the
three starter errors remain in `src/Compute.jl`.

The reference code, sample answers, and generated reference outputs are kept
outside the student archive.
