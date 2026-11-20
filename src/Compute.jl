# Student implementation | Complete the three functions in Parts 1 and 2.
# Already written for you: Network.jl reads a network file with load_network(path) and holds
# state_matrices, bound_fraction, and inducer; Reporting.jl runs the analyses and writes the
# plots and tables.

# Do not change the function names, arguments, or return types.
# Replace each placeholder error(...) with your implementation and return value.
# The state x lists the mRNA levels of the N genes, then their protein levels, in nmol/gDW.
# For the memory motif, x = [m_X, m_Y, m_Z, p_X, p_Y, p_Z]. The input u lists the control of
# each gene's transcription, a number from 0 to 1. The inducer I is in µM, and time is in min.

# ----------------------------------------------------------------------------
# Part 1a: The control of one gene
# ----------------------------------------------------------------------------

"""
    control(network::MyGeneNetwork, j::Int64, x::Array{Float64,1}, I::Float64) -> Float64

Compute the control ū_j of gene `j`: the fraction of the gene's full transcription rate,
`Vmax_X`, at the state `x` and the inducer concentration `I`.

The control comes from the promoter states of gene `j`, `network.genes[j].promoter`. The weight
of state s is its `W` times the bound fraction of each regulator `b` bound in that state:

    weight_s = W_s × f_1 × f_2 × ...,   where f = bound_fraction(network, b, x, I)

A state with no bound regulators has the weight `W_s`. The bare promoter, with nothing bound,
has the weight 1 and is not listed in the file. The control is the summed weight of the active
states divided by 1 plus the summed weight of all the listed states:

    ū_j = (sum of weight_s over the active states) / (1 + sum of weight_s over all states)

A gene with no promoter states is always transcribed at its full rate, so its control is `1.0`.

### Arguments
- `network::MyGeneNetwork`: the network, from `load_network(path)`.
- `j::Int64`: the gene, from 1 to the number of genes `N`.
- `x::Array{Float64,1}`: the state, of length `2N`: the mRNA levels, then the protein levels
  (nmol/gDW).
- `I::Float64`: the inducer concentration (µM).

### Returns
The control, a `Float64` from 0 to 1. Leave `x` unchanged.

### Errors
Throws an `ArgumentError` if `j` is not from 1 to `N`, or if `x` does not have length `2N`.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
x = zeros(6);                # no mRNA and no protein
control(network, 1, x, 0.0)  # 1.0: gene X has no promoter states
control(network, 2, x, 0.0)  # 0.005/1.005: only the leak state of gene Y has weight
```
"""
function control(network::MyGeneNetwork, j::Int64, x::Array{Float64,1}, I::Float64)::Float64

    error("Complete control in src/Compute.jl");
end

# ----------------------------------------------------------------------------
# Part 1b: The bilinear discretization
# ----------------------------------------------------------------------------

"""
    discretize(A::Array{Float64,2}, B::Array{Float64,2}, Δt::Float64) -> Tuple{Array{Float64,2},Array{Float64,2}}

Compute the bilinear discretization `(Ā, B̄)` of the state space model `dx/dt = A x + B u` for the
time step `Δt` (min):

    Ā = (I - (Δt/2) A)⁻¹ (I + (Δt/2) A)
    B̄ = (I - (Δt/2) A)⁻¹ Δt B

In this function, `I` is the identity matrix (`I` from LinearAlgebra), not the inducer. The
recursion `x_k = Ā x_(k-1) + B̄ u_k` then steps the state forward by `Δt`. Solve with the
backslash operator instead of forming the inverse of `I - (Δt/2) A`.

### Arguments
- `A::Array{Float64,2}`: the square state matrix, `n × n`.
- `B::Array{Float64,2}`: the input matrix, `n × m`.
- `Δt::Float64`: the time step (min), greater than 0.

### Returns
The tuple `(Ā, B̄)` of new matrices: `Ā` is `n × n` and `B̄` is `n × m`. Leave `A` and `B`
unchanged.

### Errors
Throws an `ArgumentError` if `Δt ≤ 0`. Throws a `DimensionMismatch` if `A` is not square or if
`B` does not have as many rows as `A`.

### Example
```julia
(A, B) = state_matrices(load_network(joinpath("data", "memory.toml")));
(Ā, B̄) = discretize(A, B, 0.1); # steps of 0.1 min
```
"""
function discretize(A::Array{Float64,2}, B::Array{Float64,2}, Δt::Float64)::Tuple{Array{Float64,2},Array{Float64,2}}

    error("Complete discretize in src/Compute.jl");
end

# ----------------------------------------------------------------------------
# Part 2: The recursion
# ----------------------------------------------------------------------------

"""
    simulate(network::MyGeneNetwork, x0::Array{Float64,1}, pulses::Array{Tuple{Float64,Float64,Float64},1};
        Δt::Float64, T::Float64) -> NamedTuple

Simulate the network from the state `x0` at time 0 to the time `T` (min) with the bilinear
recursion and the time step `Δt` (min), while the inducer follows `pulses`.

Carry out these steps, in this order:

1. Check the arguments (see Errors) before any other work.
2. Build `A` and `B` with `state_matrices(network)`, then `(Ā, B̄)` with your
   `discretize(A, B, Δt)`.
3. The number of steps is `K = round(Int64, T / Δt)`, and the times are
   `t = [0, Δt, 2Δt, ..., KΔt]`, with `t[k + 1] = k * Δt`.
4. Store `x0` as column 1 of the states. Then, for each step `k = 1, 2, ..., K`:
   - the inducer over the step is its value at the start of the step,
     `I_k = inducer(pulses, t[k])`;
   - the control of each gene `j = 1, 2, ..., N` is `u[j] = control(network, j, x, I_k)`, where
     `x` is the state at the start of the step (column `k`);
   - the state at the end of the step (column `k + 1`) is `Ā x + B̄ u`.

The control is held over each step: it uses the state at the start of the step.

### Arguments
- `network::MyGeneNetwork`: the network, from `load_network(path)`.
- `x0::Array{Float64,1}`: the state at time 0, of length `2N` (nmol/gDW).
- `pulses::Array{Tuple{Float64,Float64,Float64},1}`: the inducer pulses, each
  `(start, stop, level)` in min, min, and µM (see `inducer`). For no inducer, pass
  `Tuple{Float64,Float64,Float64}[]`.
- `Δt::Float64`: the time step (min). It has no default.
- `T::Float64`: the final time (min), a whole number of steps. It has no default.

### Returns
A named tuple `(t = t, x = X)`:
- `t::Array{Float64,1}`: the `K + 1` times (min), from `0.0` to `T`.
- `x::Array{Float64,2}`: the states, `2N × (K + 1)`. Column `k` is the state at time `t[k]`, and
  column 1 is `x0`.
Leave `x0` and `pulses` unchanged.

### Errors
Throws an `ArgumentError`, before any other work, if `Δt ≤ 0`, if `T < Δt`, if `T` is not a
whole number of steps (`abs(T / Δt - K) > 1e-9`), or if `x0` does not have length `2N`.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
x0 = zeros(6);                # start with no mRNA and no protein
pulses = [(60.0, 90.0, 5.0)]; # 5 µM of inducer from 60 to 90 min
result = simulate(network, x0, pulses; Δt = 0.1, T = 900.0);
result.x[5, end]              # protein Y at 900 min
```
"""
function simulate(network::MyGeneNetwork, x0::Array{Float64,1}, pulses::Array{Tuple{Float64,Float64,Float64},1};
    Δt::Float64, T::Float64)::NamedTuple

    error("Complete simulate in src/Compute.jl");
end
