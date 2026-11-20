# Problem Set 6 (PS6): A gene circuit that remembers a signal

A biotechnology company builds a cell-based sensor that records whether a sample
ever contained a chemical signal. The sensor cells are bacteria that carry a small
gene circuit. When the signal molecule enters a cell, the circuit switches on two
genes, Y and Z, and keeps them on after the signal is gone. Protein Z carries a
fluorescent tag, so the company reads the record as the fluorescence of the cells
hours later.

Problem Set 6 (PS6) is due: __11:59 PM ET, Saturday, December 5, 2026__. Released November 21, 2026.

## The question

The company asks:

> **How long must the signal last for the cells to remember it, and why do the
> cells remember it at all?**

The circuit is the **regulated feedback** motif, a memory motif found in the gene
networks that control cell development (Alon, 2019; see [Sources](#sources)). It
has three genes, X, Y, and Z, and one **inducer**, I: the signal molecule.

```
         I
         |  (binds protein X)
         v
         X
       /   \
      v     v
      Y <-> Z
```

- Each gene is transcribed into its mRNA, and each mRNA is translated into its
  protein. The mRNA and the protein are lost by degradation and diluted as the
  cells grow and divide.
- Proteins X, Y, and Z are **transcription factors**: proteins that bind the
  **promoter** of a gene, the DNA region where RNA polymerase starts
  transcription, and change how often the gene is transcribed.
- Gene X is always transcribed. Only protein X that carries the inducer can bind
  a promoter. Bound protein X activates genes Y and Z.
- Protein Y activates gene Z, and protein Z activates gene Y. These two
  activations form the **Y–Z loop**.

We'll model the circuit as a state space model with six states: the mRNA and
protein levels of genes X, Y, and Z. The input is the **control** of each gene, a
number from `0` to `1` that gives the fraction of the gene's full transcription
rate. We'll simulate the model with the bilinear discretization, which turns it
into a recursion that steps the state forward by a fixed time step.

You will write the control of a gene, the bilinear discretization, and the
recursion. You will then use them to find how long a pulse of inducer must last
for the cells to remember it.

## Setup

Use Julia `1.12.7`. Open this folder in VS Code, open a terminal, and run:

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

The first command installs the packages. The second runs the 40 tests, which
fail until you finish Parts 1 and 2.

Write your code in [src/Compute.jl](src/Compute.jl) and your answers in
[responses.md](responses.md). Do not edit any other files, with one exception: if
you put helper functions in a new file in `src/`, add that file to
[Include.jl](Include.jl), so that `Compute.jl` can use the helper functions.

## The network file

The circuit and its parameters are in [data/memory.toml](data/memory.toml). The
file [data/memory-no-loop.toml](data/memory-no-loop.toml) holds the same circuit
with the Y–Z loop cut: proteins Y and Z no longer activate each other's genes.
Calling `network = load_network(path)` ([src/Network.jl](src/Network.jl)) reads a
network file into `network`.

Levels inside the cells are in nmol/gDW, nanomoles per gram of dry cell weight.
The inducer concentration is in µM, in the medium around the cells. Times are in
minutes. The [data notes](data/README.md) give the sources of the values.

| Field | Meaning |
|---|---|
| `mRNA_half_life`, `protein_half_life` | The half-life of each mRNA and each protein to degradation (min). |
| `doubling_time` | The time the cells take to double (min), which sets the dilution by growth. |
| `Vmax_X` | The transcription rate of a gene when its control is `1` (nmol/gDW/min). |
| `k_L` | The translation rate constant: the proteins made per mRNA per minute (1/min). |
| `W`, `active` | The weight of a promoter state, and whether RNA polymerase transcribes the gene in that state. |
| `regulator`, `K`, `n` | The transcription factor bound in a promoter state, the level at which half of the promoters have it bound (nmol/gDW), and the binding order. |
| `inducer`, `K_I`, `n_I` | The inducer the transcription factor must carry to bind, the inducer concentration at which half of the transcription factor carries it (µM), and the binding order of the inducer. |

A **promoter state** is one arrangement of proteins bound to the promoter. The
promoter of gene Y has three listed states: RNA polymerase bound alone, which
gives a small leak of transcription; protein X bound; and protein Z bound. The
promoter of gene Z has the same states, with protein Y in place of protein Z.
Gene X has no promoter states, so its control is always `1`.

The state of the circuit is the vector `x = [m_X, m_Y, m_Z, p_X, p_Y, p_Z]`: the
three mRNA levels, then the three protein levels. The supplied code also gives you
these functions:

- `state_matrices(network)` returns the matrices `A` and `B` of the state space
  model.
- `bound_fraction(network, b, x, I)` returns the fraction of promoters that have
  the regulator `b` bound.
- `inducer(pulses, t)` returns the inducer concentration at time `t`. Each pulse
  is a tuple `(start, stop, level)`: for example, `(60.0, 90.0, 5.0)` is 5 µM of
  inducer from 60 to 90 min.

## Part 1: The control and the discretization

Implement in [src/Compute.jl](src/Compute.jl):

1. **`control(network, j, x, I)`:** returns the control of gene `j` at the state
   `x` and the inducer concentration `I`, from the weights of the gene's promoter
   states.
2. **`discretize(A, B, Δt)`:** returns the matrices `(Ā, B̄)` of the bilinear
   discretization of the state space model for the time step `Δt`.

The description above each function in [src/Compute.jl](src/Compute.jl) gives the
formula, the inputs it must reject, and the error to throw. Test with:

```bash
julia --startup-file=no --project=. testme_part_1.jl
```

## Part 2: The recursion

Part 2 uses your Part 1 functions, so finish Part 1 first.

Implement in [src/Compute.jl](src/Compute.jl):

1. **`simulate(network, x0, pulses; Δt, T)`:** runs the recursion from the state
   `x0` at time `0` to the time `T`, while the inducer follows `pulses`, and
   returns the times and the states. Each step uses the inducer concentration and
   the control from the start of the step.

When `simulate` is done, run:

```bash
julia --startup-file=no --project=. testme_part_2.jl
julia --startup-file=no --project=. runmotif.jl
```

## Run the analyses

With no inducer, the memory motif has three **steady states**, where every level
stays constant: a low state, a middle state, and a high state. The low and high
states are stable. The middle state is not: protein Y moves away from it, toward
the low state or toward the high state. Every run in `runmotif.jl` starts at the
low state of the memory motif at 0 min and ends at 900 min, with a time step of
0.1 min. Each pulse of inducer starts at 60 min. Genes Y and Z have the same
parameters and start at the same levels, so protein Y also gives the level of
protein Z, which the company reads as fluorescence. A pulse is **remembered** if
protein Y is above the middle state at 900 min, and **forgotten** otherwise.

The script `runmotif.jl` uses your three functions to run the analyses, prints the
main results, and writes these files to the `outputs/` directory:

| File | Contents |
|---|---|
| `steady-states.csv` | The protein levels at the low, middle, and high states, and the largest real part of the eigenvalues at each state. A negative value means the state is stable. |
| `memory.csv` | The protein levels every minute from 0 to 900 min in three runs: a 60-min pulse of 5 µM in the memory motif, the same pulse in the circuit with the Y–Z loop cut, and a 10-min pulse of 5 µM in the memory motif. |
| `memory.png` | Protein Y in the three runs of `memory.csv`, with the 60-min pulse shaded and a dashed line at the middle state. |
| `pulse-threshold.csv` | For each inducer level (2, 5, 10, 20, and 50 µM): the shortest remembered pulse, protein Y at the end of that pulse, the minutes from the end of that pulse until protein Y crosses the middle state, the longest forgotten pulse, and the peak of protein Y after the longest forgotten pulse. |

Use these files to answer the two questions in [responses.md](responses.md).

## Check and submit

Run the checker:

```bash
julia --startup-file=no --project=. check_submission.jl
```

The checker runs the 40 tests, checks that you answered both questions, and
writes `MANIFEST.txt`, a record of your results and files. It also deletes the
old reports in `outputs/` and, if every test passes, makes new ones as
`runmotif.jl` does. It does not submit anything.

Zip the whole folder, including `MANIFEST.txt` and `outputs/`, as
`CHEME-4800-5800-PS6-<netid>.zip` (for example,
`CHEME-4800-5800-PS6-abc123.zip`) and upload it to Canvas. If some tests still
fail at the deadline, submit anyway for partial credit. See
[RUBRIC.md](RUBRIC.md) for the grading scale.

## Policies

- **Revisions:** Submit something by the deadline, even if it does not run. You
  may then revise and resubmit until the end of the semester; your highest score
  counts. No initial submission means a score of `0` and no revisions.
- **Reference solution:** Released after the deadline. Use it to understand and
  debug your work; copying it results in a score of `0`.
- **Collaboration and AI:** Discuss ideas with classmates, but do not share code
  or solutions. You may use Julia documentation, AI tools, and internet
  resources. You must be able to explain your code.

## Sources

The regulated feedback motif comes from **Alon (2019), An Introduction to Systems
Biology: Design Principles of Biological Circuits, 2nd edition, Chapman and
Hall/CRC, chapter 5, "Positive Feedback, Bistability and Memory"**
([DOI: 10.1201/9780429283321](https://doi.org/10.1201/9780429283321)), and
**Alon (2007), "Network motifs: theory and experimental approaches," Nature
Reviews Genetics, 8(6), 450–461**
([DOI: 10.1038/nrg2102](https://doi.org/10.1038/nrg2102)). The bilinear
discretization of a state space model is used in **Gu, Goel, and Ré (2022),
"Efficiently Modeling Long Sequences with Structured State Spaces," International
Conference on Learning Representations (ICLR)**
([arXiv: 2111.00396](https://arxiv.org/abs/2111.00396)).

The [data notes](data/README.md) give the sources of the rate constants and
half-lives. The company and its sensor are made up for this problem set.
