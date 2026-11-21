# ----------------------------------------------------------------------------
# Supplied network code | Read a gene network from a TOML file.
# ----------------------------------------------------------------------------

# The gas constant R (kJ/mol/K), for the weight W = exp(−dG/RT) of each promoter state -
const GAS_CONSTANT = 8.314462618e-3;

"""
    MyBoundRegulator

One regulator bound in a promoter state.

### Fields
- `regulator::String`: the name of the gene whose protein binds, such as `"X"`.
- `K::Float64`: the available regulator level at which the bound fraction is one half
  (nmol/gDW).
- `n::Float64`: the binding order (dimensionless).
- `inducer::Bool`: `true` when the inducer I sets the fraction of the regulator that can bind.
- `K_I::Float64`: the inducer concentration at which half of the regulator carries the inducer
  (µM). `0.0` when `inducer` is `false`.
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
- `bound::Array{MyBoundRegulator,1}`: the regulators bound in this state. A state with no bound
  regulators, such as RNA polymerase bound alone, has an empty list.
"""
struct MyPromoterState
    dG::Float64
    W::Float64
    active::Bool
    bound::Array{MyBoundRegulator,1}
end

"""
    MyGene

One gene of the network.

### Fields
- `name::String`: the gene name, such as `"Y"`.
- `G::Float64`: the gene level (nmol/gDW).
- `Vmax_X::Float64`: the transcription kinetic limit (nmol/gDW/min).
- `K_X::Float64`: the transcription saturation constant (nmol/gDW).
- `Vmax_L::Float64`: the translation kinetic limit (nmol/gDW/min).
- `K_L::Float64`: the translation saturation constant (nmol/gDW).
- `promoter::Array{MyPromoterState,1}`: the promoter states, without the empty promoter.
"""
struct MyGene
    name::String
    G::Float64
    Vmax_X::Float64
    K_X::Float64
    Vmax_L::Float64
    K_L::Float64
    promoter::Array{MyPromoterState,1}
end

"""
    MyGeneNetwork

A gene network read from a TOML file by `load_network`.

### Fields
- `description::String`: one line from the file's `[metadata]` table (`""` if there is none).
- `genes::Array{MyGene,1}`: the genes, in file order. Gene `j` has mRNA state `j` and protein
  state `N + j`, where `N` is the number of genes.
- `index::Dict{String,Int64}`: the position of each gene, by name.
- `temperature::Float64`: the absolute temperature (K).
- `θ_m::Float64`, `θ_p::Float64`: the degradation rate constants of mRNA and protein (1/min).
- `μ::Float64`: the growth rate (1/min). Growth dilutes mRNA and protein at this rate.
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
    rate_from_half_life(time::Float64, location::String) -> Float64

Return the first-order rate constant `ln 2 / time` (1/min) for a half-life or doubling time in
minutes. `location` is the name of the field in the file, for the error message. Throws an
`ArgumentError` if `time` is not positive.
"""
function rate_from_half_life(time::Float64, location::String)::Float64

    # Check the time -
    if time <= 0.0
        throw(ArgumentError("$location must be positive, got $time"));
    end

    # The rate constant -
    return log(2) / time;
end

"""
    load_network(path::String) -> MyGeneNetwork

Read a gene network from the TOML file at `path`, such as `data/memory.toml`.

The file has a `[parameters]` table with `mRNA_half_life`, `protein_half_life`, and
`doubling_time` (min), and `temperature` (K). Each `[[genes]]` entry has `name`, `G`, `Vmax_X`,
`K_X`, `Vmax_L`, and `K_L`, and lists its promoter states as `[[genes.promoter]]` entries with
`dG` (kJ/mol) and `active`. The loader stores the weight `W = exp(−dG/RT)` of each state. Each
state lists its bound regulators as `[[genes.promoter.bound]]` entries with `regulator`, `K`, and
`n`. A bound regulator that carries the inducer also has `inducer = "I"`,
`inducer_effect = "activates"`, `K_I`, and `n_I`, and the file then needs an `[inputs.I]` table
with `initial = 0.0`, because there is no inducer before the first pulse.

### Errors
Throws an `ArgumentError` that names the field if a field is missing, if a half-life, the doubling
time, or the temperature is not positive, if a gene name appears twice, if a state names a
regulator that is not a gene of the network, if a bound regulator names an inducer other than
`"I"` or an effect other than `"activates"`, or if the file has no `[inputs.I]` table with
`initial = 0.0` when a regulator carries the inducer.

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

    # The loss rate constants (1/min) from the half-lives and the doubling time (min):
    # degradation of mRNA and protein, and dilution by growth -
    parameters = required_value(data, "parameters", "the file");
    mRNA_half_life = Float64(required_value(parameters, "mRNA_half_life", "[parameters]"));
    θ_m = rate_from_half_life(mRNA_half_life, "mRNA_half_life");
    protein_half_life = Float64(required_value(parameters, "protein_half_life", "[parameters]"));
    θ_p = rate_from_half_life(protein_half_life, "protein_half_life");
    doubling_time = Float64(required_value(parameters, "doubling_time", "[parameters]"));
    μ = rate_from_half_life(doubling_time, "doubling_time");

    # The temperature (K) and the thermal energy RT (kJ/mol), for the weights of the states -
    temperature = Float64(required_value(parameters, "temperature", "[parameters]"));
    if isfinite(temperature) == false || temperature <= 0.0
        throw(ArgumentError("temperature in [parameters] must be a positive number (K), got $temperature"));
    end
    RT = GAS_CONSTANT * temperature;

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
        G = Float64(required_value(table, "G", location));
        Vmax_X = Float64(required_value(table, "Vmax_X", location));
        K_X = Float64(required_value(table, "K_X", location));
        Vmax_L = Float64(required_value(table, "Vmax_L", location));
        K_L = Float64(required_value(table, "K_L", location));

        # The promoter states, each with its weight W = exp(−dG/RT) -
        states = MyPromoterState[];
        state_tables = get(table, "promoter", Dict{String,Any}[]);
        for s in 1:length(state_tables)
            location = "promoter state $s of gene \"$name\"";
            dG = Float64(required_value(state_tables[s], "dG", location));
            W = exp(-dG / RT);
            active = Bool(required_value(state_tables[s], "active", location));

            # The regulators bound in this state -
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

                # The inducer, if this regulator carries it -
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
        push!(genes, MyGene(name, G, Vmax_X, K_X, Vmax_L, K_L, states));
    end

    # A regulator that carries the inducer needs the [inputs.I] table, with no inducer at the start -
    if uses_inducer == true
        if haskey(data, "inputs") == false || haskey(data["inputs"], "I") == false
            throw(ArgumentError("a regulator carries the inducer \"I\", so the file needs an [inputs.I] table with initial = 0.0"));
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

Return the bound fraction `aⁿ / (Kⁿ + aⁿ)` of the regulator `b` at the state `x` and the inducer
concentration `I` (µM). The bound fraction multiplies the weight of each promoter state that
lists `b`.

The available level `a` is the protein level of the regulator, read from `x`. When the regulator
carries the inducer (`b.inducer` is `true`), only the fraction of the protein that carries the
inducer can bind, so `a` is the protein level times `I^n_I / (K_I^n_I + I^n_I)`. A negative level
binds nothing: the function returns `0.0`.

### Errors
Throws an `ArgumentError` if `I` is negative.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
b = network.genes[2].promoter[3].bound[1]; # protein Z bound to the promoter of Y
x = zeros(6);
x[6] = b.K;                                 # protein Z at the level K
f = bound_fraction(network, b, x, 0.0);     # 0.5
```
"""
function bound_fraction(network::MyGeneNetwork, b::MyBoundRegulator, x::Array{Float64,1}, I::Float64)::Float64

    # Check the inducer concentration -
    if I < 0.0
        throw(ArgumentError("the inducer concentration must not be negative, got $I"));
    end

    # The protein level of the regulator: protein states follow the N mRNA states -
    N = length(network.genes);
    a = x[N + network.index[b.regulator]];

    # Only the protein that carries the inducer can bind -
    if b.inducer == true
        a = a * (I^b.n_I / (b.K_I^b.n_I + I^b.n_I));
    end

    # A negative level binds nothing -
    if a <= 0.0
        return 0.0;
    end

    # The bound fraction -
    return a^b.n / (b.K^b.n + a^b.n);
end

"""
    inducer(pulses::Array{Tuple{Float64,Float64,Float64},1}, t::Float64) -> Float64

Return the inducer concentration I (µM) at time `t` (min). Each pulse is a tuple
`(start, stop, level)`: the inducer is at `level` from `start` up to, but not including, `stop`.
Overlapping pulses add. With no pulses (`Tuple{Float64,Float64,Float64}[]`), the inducer is `0.0`.

### Errors
Throws an `ArgumentError` if a pulse stops at or before it starts, or has a negative level.

### Example
```julia
pulses = [(60.0, 90.0, 5.0)]; # 5 µM from 60 to 90 min
inducer(pulses, 60.0) # 5.0
inducer(pulses, 90.0) # 0.0
```
"""
function inducer(pulses::Array{Tuple{Float64,Float64,Float64},1}, t::Float64)::Float64

    # Add the levels of the pulses that are on at time t -
    I = 0.0;
    for k in 1:length(pulses)
        (start, stop, level) = pulses[k];
        if stop <= start
            throw(ArgumentError("pulse $k stops at $stop min, which is not after its start at $start min"));
        end
        if level < 0.0
            throw(ArgumentError("pulse $k has a negative level, $level µM"));
        end
        if start <= t && t < stop
            I = I + level;
        end
    end

    # Return the concentration -
    return I;
end
