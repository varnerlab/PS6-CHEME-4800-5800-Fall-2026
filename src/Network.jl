# ----------------------------------------------------------------------------
# Supplied network code | Read a gene network from a TOML file and build its model.
# ----------------------------------------------------------------------------

"""
    MyBoundRegulator

One regulator bound in a promoter state.

### Fields
- `regulator::String`: the name of the gene whose protein binds, such as `"X"`.
- `K::Float64`: the level of the available regulator at which half of the promoters have it
  bound (nmol/gDW).
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
- `W::Float64`: the weight of the state when its regulators are fully bound (dimensionless).
- `active::Bool`: `true` when RNA polymerase transcribes the gene in this state.
- `bound::Array{MyBoundRegulator,1}`: the regulators bound in this state. A state with no bound
  regulators, such as RNA polymerase bound alone, has an empty list.
"""
struct MyPromoterState
    W::Float64
    active::Bool
    bound::Array{MyBoundRegulator,1}
end

"""
    MyGene

One gene of the network.

### Fields
- `name::String`: the gene name, such as `"Y"`.
- `Vmax_X::Float64`: the transcription rate when the gene is fully on (nmol/gDW/min).
- `k_L::Float64`: the translation rate constant: proteins made per mRNA per minute (1/min).
- `promoter::Array{MyPromoterState,1}`: the promoter states. A gene with no promoter states is
  always transcribed at `Vmax_X`.
"""
struct MyGene
    name::String
    Vmax_X::Float64
    k_L::Float64
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
- `θ_m::Float64`, `θ_p::Float64`: the degradation rate constants of mRNA and protein (1/min).
- `μ::Float64`: the growth rate (1/min). Growth dilutes mRNA and protein at this rate.
"""
struct MyGeneNetwork
    description::String
    genes::Array{MyGene,1}
    index::Dict{String,Int64}
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
`doubling_time` (min), and one `[[genes]]` entry per gene with `name`, `Vmax_X`, and `k_L`. Each
gene lists its promoter states as `[[genes.promoter]]` entries with `W` and `active`, and each
state lists its bound regulators as `[[genes.promoter.bound]]` entries with `regulator`, `K`, and
`n`. A bound regulator that carries the inducer also has `inducer = "I"`,
`inducer_effect = "activates"`, `K_I`, and `n_I`.

### Errors
Throws an `ArgumentError` that names the field if a field is missing, if a half-life or the
doubling time is not positive, if a gene name appears twice, if a state names a regulator that
is not a gene of the network, or if a bound regulator names an inducer other than `"I"` or an
effect other than `"activates"`.

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
    for j in 1:length(gene_tables)
        table = gene_tables[j];
        name = string(table["name"]);
        Vmax_X = Float64(required_value(table, "Vmax_X", "gene \"$name\""));
        k_L = Float64(required_value(table, "k_L", "gene \"$name\""));

        # The promoter states; a gene without any is always transcribed at Vmax_X -
        states = MyPromoterState[];
        state_tables = get(table, "promoter", Dict{String,Any}[]);
        for s in 1:length(state_tables)
            location = "promoter state $s of gene \"$name\"";
            W = Float64(required_value(state_tables[s], "W", location));
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
                end
                push!(bound, MyBoundRegulator(regulator, K, n, has_inducer, K_I, n_I));
            end
            push!(states, MyPromoterState(W, active, bound));
        end
        push!(genes, MyGene(name, Vmax_X, k_L, states));
    end

    # Return the network -
    return MyGeneNetwork(description, genes, index, θ_m, θ_p, μ);
end

"""
    state_matrices(network::MyGeneNetwork) -> Tuple{Array{Float64,2},Array{Float64,2}}

Return the matrices `A` and `B` of the state space model `dx/dt = A x + B u`, where the state `x`
lists the mRNA levels of the `N` genes, then their protein levels, and the input `u` lists the
control of each gene's transcription.

### Returns
- `A::Array{Float64,2}` (`2N × 2N`): each mRNA is lost at the rate `θ_m + μ`, each protein is made
  from its mRNA at the rate `k_L` and lost at the rate `θ_p + μ`.
- `B::Array{Float64,2}` (`2N × N`): gene `j` is transcribed at `Vmax_X` times its control, so
  `B[j, j] = Vmax_X` of gene `j`; every other entry is zero.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
(A, B) = state_matrices(network); # A is 6 × 6 and B is 6 × 3
```
"""
function state_matrices(network::MyGeneNetwork)::Tuple{Array{Float64,2},Array{Float64,2}}

    # Sizes -
    N = length(network.genes);
    A = zeros(2 * N, 2 * N);
    B = zeros(2 * N, N);

    # Fill in each gene's mRNA and protein rows -
    for j in 1:N
        gene = network.genes[j];
        A[j, j] = -(network.θ_m + network.μ);     # mRNA lost to degradation and dilution
        A[N + j, j] = gene.k_L;                     # protein made from its mRNA
        A[N + j, N + j] = -(network.θ_p + network.μ); # protein lost to degradation and dilution
        B[j, j] = gene.Vmax_X;                      # transcription, scaled by the control
    end

    # Return the pair -
    return (A, B);
end

"""
    bound_fraction(network::MyGeneNetwork, b::MyBoundRegulator, x::Array{Float64,1}, I::Float64) -> Float64

Return the fraction of promoters that have the regulator `b` bound, `aⁿ / (Kⁿ + aⁿ)`, at the state
`x` and the inducer concentration `I` (µM).

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
x = [0.0, 0.0, 0.0, 0.0, 0.0, 3.0];         # protein Z at its half-binding level
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
