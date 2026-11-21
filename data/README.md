# Network files

| File | Contents |
|---|---|
| [memory.toml](memory.toml) | The regulated feedback (memory) motif: the inducer I activates protein X, protein X activates genes Y and Z, and proteins Y and Z activate each other's genes. |
| [memory-no-loop.toml](memory-no-loop.toml) | The same circuit with the Y–Z loop cut: the promoter states in which protein Y or protein Z is bound are left out. Every other value matches `memory.toml`. |

The README section [The network file](../README.md#the-network-file) describes
each field. The comments in each file give the meaning and units of each number.

## Where the values come from

The values describe a bacterium such as *Escherichia coli* growing at 37 °C.
Levels inside the cells are per gram of dry cell weight (gDW). One molecule per
cell is 5.93 × 10⁻³ nmol/gDW, for a dry mass of 2.8 × 10⁻¹³ g per cell
(BioNumbers ID 103904; Milo et al., 2010, *Nucleic Acids Research*, 38,
D750–D753, [DOI: 10.1093/nar/gkp889](https://doi.org/10.1093/nar/gkp889)).

| Value | Where it comes from |
|---|---|
| `G` = 0.005 nmol/gDW, `Vmax_X` = 0.04 nmol/gDW/min, `K_X` = 0.005 nmol/gDW, `Vmax_L` = 0.05 nmol/gDW/min, `K_L` = 0.05 nmol/gDW | Values from the course materials. `G` is about one gene copy per cell. |
| mRNA half-life 2.5 min, protein half-life 35 min, doubling time 60 min | Values from the course materials for a bacterium growing in culture. |
| `temperature` = 310.15 K | 37 °C. |
| `dG` = 13.66 kJ/mol for the leak states, `dG` = −5.94 kJ/mol for the states with protein Y or protein Z bound | Values from the course materials for a leaky promoter and an activated promoter. |
| `dG` = −5.94 kJ/mol for the one promoter state of gene X | The weight of an activated promoter, so that gene X is always transcribed. |
| `dG`, `K`, and `n` of the states with protein X bound; `K` and `n` of the states with protein Y or protein Z bound | Chosen so that the circuit has a low, a middle, and a high steady state with no inducer, and so that a pulse of inducer can switch the circuit to the high state. |
| `K_I` = 10 µM, `n_I` = 1 | Values from the course materials: half of protein X carries the inducer at 10 µM. |

At the high steady state, genes Y and Z each have about 10 mRNA and about 150
protein molecules per cell.
