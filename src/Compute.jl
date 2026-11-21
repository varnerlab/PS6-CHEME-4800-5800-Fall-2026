# Student implementation | Complete the three functions in Parts 1 and 2.
# Already written for you: Network.jl reads a network file with load_network(path) and holds
# bound_fraction and inducer; Reporting.jl runs the analyses and writes the plots and tables.

# Do not change the function names, arguments, or return types.
# Replace each placeholder error(...) with your implementation and return value.
# The state x lists the mRNA levels of the N genes, then their protein levels, in nmol/gDW.
# For the memory motif, x = [m_X, m_Y, m_Z, p_X, p_Y, p_Z]. The inducer I is in µM, and time is
# in min.

# ----------------------------------------------------------------------------
# Part 1a: The control of one gene
# ----------------------------------------------------------------------------

"""
    control(network::MyGeneNetwork, j::Int64, x::Array{Float64,1}, I::Float64) -> Float64

Compute the control ū_j of gene `j`: the probability that the promoter of gene `j` is in a state
that transcribes, at the state `x` and the inducer concentration `I`.

The control comes from the promoter states of gene `j`, `network.genes[j].promoter`. The weight
of state s is its `W` times the bound fraction of each regulator `b` bound in that state:

    weight_s = W_s × f_1 × f_2 × ...,   where f = bound_fraction(network, b, x, I)

The loader stores `W_s = exp(−dG_s/RT)` for each state. A state with no bound regulators, such as
RNA polymerase bound alone, has the weight `W_s`. The empty promoter has the weight 1 and is not
listed in the file. The control is the summed weight of the active states divided by 1 plus the
summed weight of all the listed states:

    ū_j = (sum of weight_s over the active states) / (1 + sum of weight_s over all states)

A gene with no promoter states is not regulated, so its control is `1.0`.

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
control(network, 1, x, 0.0)  # about 0.909: gene X has one active state, with W about 10
control(network, 2, x, 0.0)  # about 0.005: only the leak state of gene Y has weight
```
"""
function control(network::MyGeneNetwork, j::Int64, x::Array{Float64,1}, I::Float64)::Float64

    error("Complete control in src/Compute.jl");
end

# ----------------------------------------------------------------------------
# Part 1b: The mRNA and protein balances
# ----------------------------------------------------------------------------

"""
    balances!(dx::Array{Float64,1}, x::Array{Float64,1}, parameters::NamedTuple, t::Float64) -> Nothing

Compute the rate of change of every mRNA and protein level at the state `x` and the time `t`
(min), and store it in `dx` (nmol/gDW/min). The `!` at the end of the name is a Julia convention:
the function changes its argument `dx`.

`parameters` is the named tuple `(network = network, pulses = pulses)`. The inducer
concentration at time `t` is `I = inducer(pulses, t)`. For each gene j, with mRNA level `m_j` and
protein level `p_j`:

    r_X,j = Vmax_X × G / (K_X + G)                transcription rate at a control of 1
    r_L,j = Vmax_L × m_j / (K_L + m_j)            translation rate
    dm_j/dt = r_X,j × control(network, j, x, I) − (θ_m + μ) × m_j
    dp_j/dt = r_L,j − (θ_p + μ) × p_j

The mRNA balance of gene j goes in `dx[j]`, and the protein balance in `dx[N + j]`. The argument
order (`dx`, `x`, parameters, `t`) is the one `ODEProblem` expects.

### Arguments
- `dx::Array{Float64,1}`: the rates of change, of length `2N`; the function overwrites it.
- `x::Array{Float64,1}`: the state, of length `2N`: the mRNA levels, then the protein levels
  (nmol/gDW).
- `parameters::NamedTuple`: `(network = network, pulses = pulses)`, with the network from
  `load_network(path)` and the inducer pulses as in `simulate`.
- `t::Float64`: the time (min).

### Returns
`nothing`. Leave `x` unchanged.

### Errors
Throws an `ArgumentError` if `x` or `dx` does not have length `2N`.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
parameters = (network = network, pulses = [(60.0, 120.0, 5.0)]);
dx = zeros(6);
balances!(dx, zeros(6), parameters, 70.0); # dx[4] = 0.0: protein X needs mRNA X first
```
"""
function balances!(dx::Array{Float64,1}, x::Array{Float64,1}, parameters::NamedTuple, t::Float64)::Nothing

    error("Complete balances! in src/Compute.jl");
end

# ----------------------------------------------------------------------------
# Part 2: The ODE solver
# ----------------------------------------------------------------------------

"""
    simulate(network::MyGeneNetwork, x0::Array{Float64,1}, pulses::Array{Tuple{Float64,Float64,Float64},1}; T::Float64) -> NamedTuple

Solve the mRNA and protein balances from the state `x0` at time 0 to the final time `T` (min),
with the inducer pulses `pulses`, and return the state once per minute.

Each pulse is a tuple `(start, stop, level)`: the inducer is at `level` (µM) from `start` up to,
but not including, `stop` (min). Overlapping pulses add, as in `inducer`. With no pulses, pass
`Tuple{Float64,Float64,Float64}[]`.

Build an `ODEProblem` from `balances!`, a copy of `x0`, the time span `(0.0, T)`, and the
parameters `(network = network, pulses = pulses)`. Solve it with `solve` and the `Tsit5()`
method, with these settings:

- `tstops`: every pulse start and stop that lies strictly between 0 and `T`. The inducer jumps at
  these times, so the solver must stop at each one.
- `saveat = 1.0`, so that the solver returns the state once per minute.
- `reltol = 1e-8` and `abstol = 1e-12`.

### Arguments
- `network::MyGeneNetwork`: the network, from `load_network(path)`.
- `x0::Array{Float64,1}`: the state at time 0, of length `2N` (nmol/gDW).
- `pulses::Array{Tuple{Float64,Float64,Float64},1}`: the inducer pulses.
- `T::Float64` (keyword): the final time, a positive whole number of minutes.

### Returns
The named tuple `(t = t, x = x)`:
- `t::Array{Float64,1}`: the times `0.0, 1.0, …, T`, of length `T + 1` (min).
- `x::Array{Float64,2}`: the state at each time, one column per time, of size `2N × (T + 1)`
  (nmol/gDW). Column 1 is `x0`.

Leave `x0` and `pulses` unchanged.

### Errors
Throws an `ArgumentError` if `T` is not a positive whole number of minutes, if `x0` does not have
length `2N`, or if a pulse stops at or before it starts or has a negative level. Throws an
`ErrorException` if the solver stops before `T`.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
run = simulate(network, zeros(6), [(60.0, 120.0, 5.0)]; T = 900.0);
run.x[5, end] # protein Y at 900 min (nmol/gDW)
```
"""
function simulate(network::MyGeneNetwork, x0::Array{Float64,1}, pulses::Array{Tuple{Float64,Float64,Float64},1};
    T::Float64)::NamedTuple

    error("Complete simulate in src/Compute.jl");
end
