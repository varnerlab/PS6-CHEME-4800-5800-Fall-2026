# ----------------------------------------------------------------------------
# Supplied network code | Read a gene network from a TOML file.
# ----------------------------------------------------------------------------

# The gas constant R (kJ/mol/K), for the weight W = exp(−dG/RT) of each promoter state -
const GAS_CONSTANT = 8.314462618e-3;

# Avogadro's number (1/mol), to turn copies per cell into nmol/gDW -
const AVOGADRO = 6.02214076e23;

"""
    MyBoundRegulator

One transcription factor bound in a promoter state.

### Fields
- `regulator::String`: the name of the gene whose protein binds, such as `"X"`.
- `K::Float64`: the concentration of the available transcription factor at which the bound
  fraction is one half (nmol/gDW).
- `n::Float64`: the binding order (dimensionless).
- `inducer::Bool`: `true` when the inducer I sets the fraction of the transcription factor that
  can bind.
- `K_I::Float64`: the inducer concentration at which half of the transcription factor is bound to
  the inducer (µM). `0.0` when `inducer` is `false`.
- `n_I::Float64`: the binding order of the inducer. `0.0` when `inducer` is `false`.
"""
struct MyBoundRegulator
    regulator::String
    K::Float64
    n::Float64
    inducer::Bool
    K_I::Float64
    n_I::Float64
end

"""
    MyPromoterState

One state of a gene's promoter.

### Fields
- `dG::Float64`: the Gibbs energy of the state relative to the empty promoter (kJ/mol).
- `W::Float64`: the weight of the state, `exp(−dG/RT)` (dimensionless).
- `active::Bool`: `true` when RNA polymerase transcribes the gene in this state.
- `bound::Array{MyBoundRegulator,1}`: the transcription factors bound in this state. A state
  with no transcription factor bound, such as RNA polymerase bound alone, has an empty list.
"""
struct MyPromoterState
    dG::Float64
    W::Float64
    active::Bool
    bound::Array{MyBoundRegulator,1}
end

"""
    MyGene

One gene of the network. The loader computes each number from the network file.

### Fields
- `name::String`: the gene name, such as `"Y"`.
- `G::Float64`: the gene concentration (nmol/gDW), from the number of copies of the gene per cell.
- `Vmax_X::Float64`: the transcription kinetic limit (nmol/gDW/min).
- `τ_X::Float64`: the transcription time constant (dimensionless).
- `K_X::Float64`: the transcription saturation constant (nmol/gDW).
- `Vmax_L::Float64`: the translation kinetic limit (nmol/gDW/min).
- `τ_L::Float64`: the translation time constant (dimensionless).
- `K_L::Float64`: the translation saturation constant (nmol/gDW).
- `promoter::Array{MyPromoterState,1}`: the promoter states, without the empty promoter.
"""
struct MyGene
    name::String
    G::Float64
    Vmax_X::Float64
    τ_X::Float64
    K_X::Float64
    Vmax_L::Float64
    τ_L::Float64
    K_L::Float64
    promoter::Array{MyPromoterState,1}
end

"""
    MyGeneNetwork

A gene network read from a TOML file by `load_network`.

### Fields
- `description::String`: one line from the file's `[metadata]` table (`""` if there is none).
- `genes::Array{MyGene,1}`: the genes, in file order. The mRNA concentration of gene `j` is
  entry `j` of the state `x`, and its protein concentration is entry `N + j`, where `N` is the
  number of genes.
- `index::Dict{String,Int64}`: the position of each gene, by name.
- `temperature::Float64`: the absolute temperature (K).
- `θ_m::Float64`, `θ_p::Float64`: the degradation rate constants of mRNA and protein (1/min).
- `μ::Float64`: the specific growth rate (1/min). Growth dilutes mRNA and protein at this rate.
"""
struct MyGeneNetwork
    description::String
    genes::Array{MyGene,1}
    index::Dict{String,Int64}
    temperature::Float64
    θ_m::Float64
    θ_p::Float64
    μ::Float64
end

"""
    required_value(table::Dict, key::String, location::String)

Return `table[key]`. Throws an `ArgumentError` that names the key and its location in the file
if the key is missing.
"""
function required_value(table::Dict, key::String, location::String)

    # Check that the key is there -
    if haskey(table, key) == false
        throw(ArgumentError("missing field `$key` in $location"));
    end

    # Return its value -
    return table[key];
end

"""
    required_positive(table::Dict, key::String, location::String) -> Float64

Return `table[key]` as a `Float64`. Throws an `ArgumentError` that names the key and its location
in the file if the key is missing, or if its value is not a finite number above zero.
"""
function required_positive(table::Dict, key::String, location::String)::Float64

    # Read the value -
    value = Float64(required_value(table, key, location));

    # Check that it is a finite number above zero -
    if isfinite(value) == false || value <= 0.0
        throw(ArgumentError("`$key` in $location must be a positive number, got $value"));
    end

    # Return the value -
    return value;
end

"""
    load_network(path::String) -> MyGeneNetwork

Read a gene network from the TOML file at `path`, such as `data/memory.toml`.

### The file
The `[parameters]` table describes the cell, RNA polymerase, and ribosomes:
- the cell: `mass_of_single_cell` (gDW/cell), `volume_of_single_cell` (L/cell), `doubling_time`
  (min), and `temperature` (K);
- RNA polymerase: `copies_of_rnap_per_cell`, `transcription_elongation_rate` (nt/s),
  `transcription_initiation_time_constant` (s), and `transcription_saturation_constant` (µM);
- ribosomes: `copies_of_ribosome_per_cell`, `translation_elongation_rate` (aa/s),
  `translation_initiation_time_constant` (s), `translation_saturation_factor`, and
  `average_polysome_number`;
- degradation: `mRNA_half_life` and `protein_half_life` (min).

Each `[[genes]]` entry has `name`, `gene_copy_number`, `length_of_gene_in_nt`, and
`length_of_protein_in_aa`, and lists its promoter states as `[[genes.promoter]]` entries with
`dG` (kJ/mol) and `active`. Each state lists its bound transcription factors as
`[[genes.promoter.bound]]` entries with `regulator`, `K`, and `n`. A transcription factor that
needs the inducer to bind also has
`inducer = "I"`, `inducer_effect = "activates"`, `K_I`, and `n_I`, and the file then needs an
`[inputs.I]` table with `initial = 0.0`, because there is no inducer before the first pulse.

### What the loader computes
For each gene, the loader computes the fields of `MyGene`:

    G = c × n_G
    Vmax_X = c × R_X × (60 × v_X) / l_G          τ_X = (v_X / l_G) × t_X       K_X = s × K_µM
    Vmax_L = K_P × c × R_L × (60 × v_L) / l_P    τ_L = (v_L / l_P) × t_L       K_L = F × K_X

with these symbols:
- c = 10⁹ / (N_A × `mass_of_single_cell`): one copy per cell, in nmol/gDW, where N_A is
  Avogadro's number;
- s = 10³ × `volume_of_single_cell` / `mass_of_single_cell`: 1 µM inside the cell, in nmol/gDW;
- n_G = `gene_copy_number`, l_G = `length_of_gene_in_nt`, l_P = `length_of_protein_in_aa`;
- R_X = `copies_of_rnap_per_cell`, v_X = `transcription_elongation_rate`,
  t_X = `transcription_initiation_time_constant`, K_µM = `transcription_saturation_constant`;
- R_L = `copies_of_ribosome_per_cell`, v_L = `translation_elongation_rate`,
  t_L = `translation_initiation_time_constant`, F = `translation_saturation_factor`,
  K_P = `average_polysome_number`.

The factor 60 turns the elongation rates from per second to per minute. The loader also stores
the weight W = exp(−dG/RT) of each promoter state, the degradation rate constants θ_m and θ_p
(ln 2 / half-life), and the specific growth rate μ (ln 2 / doubling time).

### Errors
Throws an `ArgumentError` that names the field if a field is missing, if a number in
`[parameters]`, a gene copy number, or a length is not positive, if a gene name appears twice, if
a state names a `regulator` that is not a gene of the network, if a transcription factor names
an inducer other than `"I"` or an effect other than `"activates"`, or if the file has no
`[inputs.I]` table with `initial = 0.0` when a transcription factor needs the inducer to bind.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
network.genes[2].name # "Y"
```
"""
function load_network(path::String)::MyGeneNetwork

    # Read the file -
    data = TOML.parsefile(path);

    # The description, if the file has one -
    description = "";
    if haskey(data, "metadata") && haskey(data["metadata"], "description")
        description = string(data["metadata"]["description"]);
    end

    # The degradation rate constants and the specific growth rate (1/min), from the half-lives
    # and the doubling time (min) -
    parameters = required_value(data, "parameters", "the file");
    θ_m = log(2) / required_positive(parameters, "mRNA_half_life", "[parameters]");
    θ_p = log(2) / required_positive(parameters, "protein_half_life", "[parameters]");
    μ = log(2) / required_positive(parameters, "doubling_time", "[parameters]");

    # The temperature (K) and the thermal energy RT (kJ/mol), for the weights of the states -
    temperature = required_positive(parameters, "temperature", "[parameters]");
    RT = GAS_CONSTANT * temperature;

    # The unit conversions from the dry mass (gDW/cell) and the volume (L/cell) of one cell:
    # one copy per cell, and 1 µM inside the cell, each in nmol/gDW -
    mass = required_positive(parameters, "mass_of_single_cell", "[parameters]");
    volume = required_positive(parameters, "volume_of_single_cell", "[parameters]");
    per_copy = 1.0e9 / (AVOGADRO * mass);
    per_µM = 1.0e3 * volume / mass;

    # RNA polymerase (nmol/gDW), its elongation rate (nt/s), the transcription initiation time
    # constant (s), and the transcription saturation constant (nmol/gDW) -
    R_X = per_copy * required_positive(parameters, "copies_of_rnap_per_cell", "[parameters]");
    v_X = required_positive(parameters, "transcription_elongation_rate", "[parameters]");
    t_X = required_positive(parameters, "transcription_initiation_time_constant", "[parameters]");
    K_µM = required_positive(parameters, "transcription_saturation_constant", "[parameters]");
    K_X = per_µM * K_µM;

    # Ribosomes (nmol/gDW), their elongation rate (aa/s), the translation initiation time constant
    # (s), the translation saturation constant (nmol/gDW), and the ribosomes on each mRNA -
    R_L = per_copy * required_positive(parameters, "copies_of_ribosome_per_cell", "[parameters]");
    v_L = required_positive(parameters, "translation_elongation_rate", "[parameters]");
    t_L = required_positive(parameters, "translation_initiation_time_constant", "[parameters]");
    K_L = required_positive(parameters, "translation_saturation_factor", "[parameters]") * K_X;
    K_P = required_positive(parameters, "average_polysome_number", "[parameters]");

    # First pass: the gene names and their positions, so that a state can name any gene -
    gene_tables = required_value(data, "genes", "the file");
    index = Dict{String,Int64}();
    for j in 1:length(gene_tables)
        name = string(required_value(gene_tables[j], "name", "gene $j"));
        if haskey(index, name)
            throw(ArgumentError("gene name \"$name\" appears twice"));
        end
        index[name] = j;
    end

    # Second pass: build each gene with its promoter states -
    genes = MyGene[];
    uses_inducer = false;
    for j in 1:length(gene_tables)
        table = gene_tables[j];
        name = string(table["name"]);
        location = "gene \"$name\"";

        # The copies of gene j per cell, and the lengths of the gene (nt) and its protein (aa) -
        copies = required_positive(table, "gene_copy_number", location);
        l_G = required_positive(table, "length_of_gene_in_nt", location);
        l_P = required_positive(table, "length_of_protein_in_aa", location);

        # The gene concentration (nmol/gDW), the kinetic limits (nmol/gDW/min), and the time
        # constants (dimensionless) of gene j. The factor 60 turns per second into per minute -
        G = per_copy * copies;
        Vmax_X = R_X * (60.0 * v_X) / l_G;
        τ_X = (v_X / l_G) * t_X;
        Vmax_L = K_P * R_L * (60.0 * v_L) / l_P;
        τ_L = (v_L / l_P) * t_L;

        # The promoter states, each with its weight W = exp(−dG/RT) -
        states = MyPromoterState[];
        state_tables = get(table, "promoter", Dict{String,Any}[]);
        for s in 1:length(state_tables)
            location = "promoter state $s of gene \"$name\"";
            dG = Float64(required_value(state_tables[s], "dG", location));
            W = exp(-dG / RT);
            active = Bool(required_value(state_tables[s], "active", location));

            # The transcription factors bound in this state -
            bound = MyBoundRegulator[];
            bound_tables = get(state_tables[s], "bound", Dict{String,Any}[]);
            for b in 1:length(bound_tables)
                entry = bound_tables[b];
                regulator = string(required_value(entry, "regulator", location));
                if haskey(index, regulator) == false
                    throw(ArgumentError("$location names regulator \"$regulator\", which is not a gene of the network"));
                end
                K = Float64(required_value(entry, "K", location));
                n = Float64(required_value(entry, "n", location));

                # The inducer, if this transcription factor needs the inducer to bind -
                has_inducer = haskey(entry, "inducer");
                K_I = 0.0;
                n_I = 0.0;
                if has_inducer == true
                    if string(entry["inducer"]) != "I"
                        throw(ArgumentError("$location names inducer \"$(entry["inducer"])\"; the only inducer is \"I\""));
                    end
                    effect = string(required_value(entry, "inducer_effect", location));
                    if effect != "activates"
                        throw(ArgumentError("$location has inducer_effect \"$effect\"; the only effect is \"activates\""));
                    end
                    K_I = Float64(required_value(entry, "K_I", location));
                    n_I = Float64(required_value(entry, "n_I", location));
                    uses_inducer = true;
                end
                push!(bound, MyBoundRegulator(regulator, K, n, has_inducer, K_I, n_I));
            end
            push!(states, MyPromoterState(dG, W, active, bound));
        end
        push!(genes, MyGene(name, G, Vmax_X, τ_X, K_X, Vmax_L, τ_L, K_L, states));
    end

    # A transcription factor that needs the inducer requires the [inputs.I] table with initial = 0.0 -
    if uses_inducer == true
        if haskey(data, "inputs") == false || haskey(data["inputs"], "I") == false
            throw(ArgumentError("a regulator needs the inducer \"I\" to bind, so the file needs an [inputs.I] table with initial = 0.0"));
        end
        initial = Float64(required_value(data["inputs"]["I"], "initial", "[inputs.I]"));
        if initial != 0.0
            throw(ArgumentError("initial in [inputs.I] must be 0.0, got $initial, because there is no inducer before the first pulse"));
        end
    end

    # Return the network -
    return MyGeneNetwork(description, genes, index, temperature, θ_m, θ_p, μ);
end

"""
    bound_fraction(network::MyGeneNetwork, b::MyBoundRegulator, x::Array{Float64,1}, I::Float64) -> Float64

Return the bound fraction `aⁿ / (Kⁿ + aⁿ)` of the transcription factor `b` at the state `x` and
the inducer concentration `I` (µM). The bound fraction multiplies the weight of each promoter state that
lists `b`.

The available concentration `a` is the protein concentration of the transcription factor, read
from `x`. When the transcription factor needs the inducer to bind (`b.inducer` is `true`), only
the fraction of the protein that is bound to the inducer can bind, so `a` is the protein concentration times
`I^n_I / (K_I^n_I + I^n_I)`. A negative concentration binds nothing: the function returns `0.0`.

### Errors
Throws an `ArgumentError` if `I` is negative.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
b = network.genes[2].promoter[3].bound[1]; # protein Z bound to the promoter of Y
x = zeros(6);
x[6] = b.K;                                 # protein Z at the concentration K
f = bound_fraction(network, b, x, 0.0);     # 0.5
```
"""
function bound_fraction(network::MyGeneNetwork, b::MyBoundRegulator, x::Array{Float64,1}, I::Float64)::Float64

    # Check the inducer concentration -
    if I < 0.0
        throw(ArgumentError("the inducer concentration must not be negative, got $I"));
    end

    # The protein concentration of the transcription factor: in x, the N protein concentrations
    # follow the N mRNA concentrations -
    N = length(network.genes);
    a = x[N + network.index[b.regulator]];

    # Only the protein that is bound to the inducer can bind -
    if b.inducer == true
        a = a * (I^b.n_I / (b.K_I^b.n_I + I^b.n_I));
    end

    # A negative concentration binds nothing -
    if a <= 0.0
        return 0.0;
    end

    # The bound fraction -
    return a^b.n / (b.K^b.n + a^b.n);
end

"""
    inducer(pulses::Array{Tuple{Float64,Float64,Float64},1}, t::Float64) -> Float64

Return the inducer concentration I (µM) at time `t` (min). Each pulse is a tuple
`(start, stop, concentration)`: the inducer concentration is `concentration` from `start` up to,
but not including, `stop`. Overlapping pulses add. With no pulses
(`Tuple{Float64,Float64,Float64}[]`), the inducer concentration is `0.0`.

### Errors
Throws an `ArgumentError` if a pulse stops at or before it starts, or has a negative concentration.

### Example
```julia
pulses = [(60.0, 90.0, 5.0)]; # 5 µM from 60 to 90 min
inducer(pulses, 60.0) # 5.0
inducer(pulses, 90.0) # 0.0
```
"""
function inducer(pulses::Array{Tuple{Float64,Float64,Float64},1}, t::Float64)::Float64

    # Add the concentrations of the pulses that are on at time t -
    I = 0.0;
    for k in 1:length(pulses)
        (start, stop, concentration) = pulses[k];
        if stop <= start
            throw(ArgumentError("pulse $k stops at $stop min, which is not after its start at $start min"));
        end
        if concentration < 0.0
            throw(ArgumentError("pulse $k has a negative concentration, $concentration µM"));
        end
        if start <= t && t < stop
            I = I + concentration;
        end
    end

    # Return the concentration -
    return I;
end
