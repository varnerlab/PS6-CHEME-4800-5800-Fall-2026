# Problem Set 6 (PS6): A gene circuit that remembers a signal

A biotechnology company builds a cell-based sensor that records whether a sample
contains a chemical signal molecule. The sensor cells are bacteria that carry a
small gene circuit. The signal molecule enters a cell and binds to the sensor
protein X, which then activates the expression of two other genes, Y and Z. Genes
Y and Z form the memory of the circuit: protein Y activates gene Z, and protein Z
activates gene Y. Protein Z carries a small fluorescent tag, so the company can
read the record of the signal by measuring the fluorescence of the cells. The
company wants to understand how long the cells must be exposed to the signal for
the cells to remember it. Let's explore the circuit and its memory.

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
- Gene X is always transcribed. Only protein X that is bound to the signal
  molecule (inducer) can bind the promoters of genes Y and Z. Bound protein X
  activates genes Y and Z.
- Protein Y activates gene Z, and protein Z activates gene Y. These two
  activations form the **Y–Z loop**.

We'll model the circuit with six material balances: one for the mRNA concentration
and one for the protein concentration of each gene. The transcription rate of each
gene is set by its kinetic limit, the gene concentration, the competition of the
other genes for RNA polymerase, and the value of its **control variable**, a
number from `0` to `1`: the probability that the gene's promoter is in an active
state, a state in which RNA polymerase transcribes the gene. The translation rate
of each mRNA rises and saturates as the mRNA concentration rises. The translation
kinetic limit and the competition of the other mRNAs for ribosomes set the
translation rate. Degradation and dilution remove mRNA and protein. We'll solve
the material balances with an ODE solver, the `Tsit5` method from the
OrdinaryDiffEq package.

You will write the control variable of a gene, the material balances, and the
function that solves the material balances. You will then use these three
functions to find how long a pulse of inducer must last for the cells to remember
the signal.

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

Intracellular concentrations are in units of nmol/gDW, nanomoles per gram of dry
cell weight. The inducer, in the medium around the cells, has concentration units
of µM, micromoles per liter. Times are in minutes. The [data notes](data/README.md)
give the sources of the values.

| Field | Meaning |
|---|---|
| `mass_of_single_cell`, `volume_of_single_cell` | The dry mass (gDW) and the volume (L) of one cell. The loader uses them to turn copies per cell and µM into nmol/gDW. |
| `doubling_time` | The time the cells take to double (min), which sets the dilution by growth. |
| `temperature` | The absolute temperature (K), which sets the weight of each promoter state. |
| `copies_of_rnap_per_cell`, `transcription_elongation_rate`, `transcription_initiation_time_constant`, `transcription_saturation_constant` | The copies of RNA polymerase per cell, the elongation rate of RNA polymerase (nt/s), the time RNA polymerase takes to start transcription (s), and the transcription saturation constant (µM). |
| `copies_of_ribosome_per_cell`, `translation_elongation_rate`, `translation_initiation_time_constant`, `translation_saturation_factor`, `average_polysome_number` | The copies of ribosomes per cell, the elongation rate of ribosomes (aa/s), the time a ribosome takes to start translation (s), the translation saturation constant divided by the transcription saturation constant, and the number of ribosomes on each mRNA. |
| `mRNA_half_life`, `protein_half_life` | The degradation half-life of each mRNA and of each protein (min). |
| `[inputs.I]`, `initial` | The inducer, and its concentration before the first pulse, `0.0` µM. |
| `gene_copy_number`, `length_of_gene_in_nt`, `length_of_protein_in_aa` | The copies of a gene per cell, the length of the gene (nt), and the length of its protein (aa). |
| `dG`, `active` | The Gibbs energy of a promoter state relative to the empty promoter (kJ/mol), and whether RNA polymerase transcribes the gene in that state. |
| `regulator`, `K`, `n` | The transcription factor bound in a promoter state, the concentration of the available transcription factor at which its bound fraction is one half (nmol/gDW), and the binding order. |
| `inducer`, `K_I`, `n_I` | The inducer the transcription factor needs to bind, the inducer concentration at which half of the transcription factor is bound to the inducer (µM), and the binding order of the inducer. |

From these values, the loader computes the gene concentration, the kinetic
limits, the saturation constants, and the time constants of each gene. The
description of `load_network` in [src/Network.jl](src/Network.jl) shows how.

A **promoter state of gene Y** is one arrangement of proteins bound to the
promoter. The promoter of gene Y has three listed states: RNA polymerase bound
alone, which gives a small leak of transcription; RNA polymerase with protein X
bound; and RNA polymerase with protein Z bound.

The **promoter of gene Z** has the same states, with protein Y in place of protein
Z. Gene X lists one state, RNA polymerase bound alone, so gene X is always
transcribed. The loader turns the Gibbs energy `dG` of each state into the state's
weight `W`: the lower the Gibbs energy of a state, the larger its weight.

The **state of the circuit** is the vector `x = [m_X, m_Y, m_Z, p_X, p_Y, p_Z]`:
the intracellular concentrations of the three mRNAs, then of the three proteins.
The supplied code also gives you these functions:

- `bound_fraction(network, b, x, I)` returns the bound fraction of `b`, a
  transcription factor bound in a promoter state, at the state `x` and the inducer
  concentration `I`. The bound fraction multiplies the weight `W` of that promoter
  state.
- `inducer(pulses, t)` returns the inducer concentration at time `t`. Each pulse
  is a tuple `(start, stop, concentration)`: the inducer concentration is
  `concentration` (µM) from `start` to `stop` (min). For example,
  `(60.0, 90.0, 5.0)` is 5 µM of inducer from 60 to 90 min.

## Part 1: The control variable and the material balances

Implement in [src/Compute.jl](src/Compute.jl):

1. **`control(network, j, x, I)`:** returns the control variable of gene `j` at
   the state `x` and the inducer concentration `I`, from the weights of the
   promoter states of gene `j`.
2. **`balances!(dx, x, parameters, t)`:** computes the rate of change of the
   intracellular concentration of every mRNA and protein at the state `x` and the
   time `t`, and stores the rates in `dx`. The argument `parameters` holds the
   network and the inducer pulses.

Read the description above each function in [src/Compute.jl](src/Compute.jl)
before you write any code. Each description explains the model behind the
function and works through an example for gene Y with numbers you can check. Each
function body has its steps written for you: write the code that each `TODO`
comment asks for, then replace the `error(...)` line at the end of the function
with the return value. Test with:

```bash
julia --startup-file=no --project=. testme_part_1.jl
```

## Part 2: The ODE solver

Part 2 uses your Part 1 functions, so finish Part 1 first.

Implement in [src/Compute.jl](src/Compute.jl):

1. **`simulate(network, x0, pulses; T)`:** solves the material balances with the
   `Tsit5` method from the initial state `x0` at time `0` to the final time `T`,
   while the inducer concentration follows `pulses`, and returns the state every
   minute. The inducer concentration jumps at the start and the stop of each
   pulse, so the solver must land exactly on those times.

The description above `simulate` explains what the ODE solver does, gives the two
lines that call the solver, explains each solver setting, and says what the solver
returns. As in Part 1, write the code that each `TODO` comment asks for.

When `simulate` is done, run:

```bash
julia --startup-file=no --project=. testme_part_2.jl
julia --startup-file=no --project=. runmotif.jl
```

## Run the analyses

With no inducer, the memory motif has three **steady states**, where every
intracellular concentration stays constant: a low state, a middle state, and a
high state. The low and high states are stable, so the memory motif is
**bistable**. The middle state is unstable: the concentration of protein Y moves
away from the middle state, toward the low state or toward the high state. Every
run in `runmotif.jl` starts at the low state of the memory motif at 0 min and ends
at 900 min, and each pulse of inducer starts at 60 min. Genes Y and Z have the
same parameters and start at the same concentrations, so the concentration of
protein Y equals the concentration of protein Z, which the company measures as
fluorescence. A pulse is **remembered** if the concentration of protein Y at 900
min is above its concentration at the middle state.

The script `runmotif.jl` uses your three functions to run the analyses, prints the
main results, and writes these files to the `outputs/` directory:

| File | Contents |
|---|---|
| `steady-states.csv` | The protein concentrations at the low, middle, and high steady states, and the largest real part of the eigenvalues of the Jacobian of the material balances at each steady state. A negative largest real part means the steady state is stable. |
| `memory.csv` | The protein concentrations every minute from 0 to 900 min in three runs: a 60-min pulse of 5 µM inducer in the memory motif, the same pulse in the circuit with the Y–Z loop cut, and a 10-min pulse of 5 µM inducer in the memory motif. |
| `memory.png` | The concentration of protein Y in the three runs of `memory.csv`, with the 60-min pulse shaded and a dashed line at the concentration of protein Y at the middle state. |
| `pulse-threshold.csv` | The shortest remembered pulse at each inducer concentration (2, 5, 10, 20, and 50 µM). The script tries pulses of 1, 2, 3, …, 600 min. |

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

Zip the whole folder, including `MANIFEST.txt` and any reports in `outputs/`, as
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
([DOI: 10.1038/nrg2102](https://doi.org/10.1038/nrg2102)).

The [data notes](data/README.md) give the sources of the parameter values. The
company and its sensor are made up for this problem set.
