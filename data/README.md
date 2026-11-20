# Network files

| File | Contents |
|---|---|
| [memory.toml](memory.toml) | The regulated feedback (memory) motif: the inducer I activates protein X, protein X activates genes Y and Z, and proteins Y and Z activate each other's genes. |
| [memory-no-loop.toml](memory-no-loop.toml) | The same circuit with the Y–Z loop cut: the promoter states in which protein Y or protein Z is bound have the weight `0`. Every other value matches `memory.toml`. |

The README section [The network file](../README.md#the-network-file) describes
each field. The comments in each file give the meaning and units of each number.

## Where the values come from

The rate constants describe a bacterium such as *Escherichia coli*. Each gene has
one copy per cell. Levels inside the cells are per gram of dry cell weight
(gDW); one molecule per cell is 5.93 × 10⁻³ nmol/gDW.

| Value | Where it comes from |
|---|---|
| `Vmax_X` = 8.475 × 10⁻³ nmol/gDW/min | One gene copy starts one transcript every 42 s, the transcription initiation time. |
| `k_L` = 12 1/min | Four ribosomes translate each mRNA at the same time, each at 16.5 amino acids per second, along a protein of 330 amino acids. |
| mRNA half-life 2.5 min, protein half-life 35 min, doubling time 60 min | Values from the course materials for a bacterium growing in culture. |
| `K`, `W`, `n` | Chosen so that the circuit has a low, a middle, and a high steady state with no inducer. Fully on, a gene makes about 5 mRNA and about 1,900 protein molecules per cell. |
| `K_I` = 10 µM, `n_I` = 1 | Chosen so that half of protein X carries the inducer at 10 µM. |

The cell constants come from the course materials, which cite BioNumbers (Milo et
al., 2010, *Nucleic Acids Research*, 38, D750–D753,
[DOI: 10.1093/nar/gkp889](https://doi.org/10.1093/nar/gkp889)):

| Constant | Value | BioNumbers ID |
|---|---|---|
| Dry mass of one cell | 2.8 × 10⁻¹³ g | 103904 |
| Transcription initiation time | 42 s | (course value) |
| Translation elongation rate | 16.5 amino acids/s | (course value) |
| Ribosomes per mRNA | 4 | (course value) |
| Protein length | 330 amino acids | (course value) |
