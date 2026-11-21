# Student implementation | Complete the three functions in Parts 1 and 2.
# Already written for you: Network.jl reads a network file with load_network(path) and holds
# bound_fraction and inducer; Reporting.jl runs the analyses and writes the plot and the tables.

# Do not change the function names, arguments, or return types.
# Read the description above each function first: it explains the model behind the function and
# works through an example. Each function body has its steps written for you. Write the code that
# each TODO comment asks for, then replace the error(...) line at the end of the function with the
# return value.
# The state x lists the mRNA concentrations of the N genes, then their protein concentrations,
# in nmol/gDW. For the memory motif, x = [m_X, m_Y, m_Z, p_X, p_Y, p_Z]. The inducer I is in µM,
# and time is in min.

# ----------------------------------------------------------------------------
# Part 1a: The control variable of one gene
# ----------------------------------------------------------------------------

"""
    control(network::MyGeneNetwork, j::Int64, x::Array{Float64,1}, I::Float64) -> Float64

Compute the control variable ū_j of gene `j` at the state `x` and the inducer concentration `I`.

### What the control variable is
RNA polymerase transcribes gene `j` only when the promoter of gene `j` is in an active state.
Transcription factors bind and leave the promoter much faster than the mRNA and protein
concentrations change, so the model uses the probability of each promoter state at the current
state `x`. The control variable ū_j is the probability that the promoter of
gene `j` is in an active state: a number from 0 to 1. In `balances!`, the transcription rate of
gene `j` is its transcription rate at a control variable of 1, times ū_j.

### The promoter states
The empty promoter, with nothing bound, is the reference state: it has the weight 1, it is not
active, and it is not listed in the network file. The listed states of gene `j` are the entries
of `network.genes[j].promoter`. Each state `s` has three fields:

- `s.W`: the weight of the state, `W = exp(−dG/RT)`. The loader computes `W` from the Gibbs
  energy `dG` of the state and the temperature: the lower the Gibbs energy, the larger the weight.
- `s.active`: `true` if RNA polymerase transcribes the gene in this state.
- `s.bound`: the list of transcription factors bound in this state, which can be empty. For each
  entry `b` of the list, `bound_fraction(network, b, x, I)` returns the bound fraction
  f_b = aⁿ / (Kⁿ + aⁿ), a number from 0 to 1. Here a is the concentration of the transcription
  factor that can bind; for protein X, only the part that is bound to the inducer can bind.

The promoter of gene Y in `data/memory.toml` lists three states, all active:

    state   bound with RNA polymerase   dG (kJ/mol)   W
    1       nothing (the leak)          13.66         0.00501
    2       protein X                   5.50          0.1185
    3       protein Z                   −5.94         10.01

### The weights and the control variable
The weight ω_s of a listed state is its `W` times the bound fraction of each transcription factor
bound in it. A state with an empty `bound` list has the weight `W`:

    ω_s = W_s × f_1 × f_2 × ...

The partition function Z_j adds the weights of all the states, with 1 for the empty promoter:

    Z_j = 1 + (the sum of ω_s over the listed states)

The probability of state s is ω_s / Z_j, so the control variable adds up the probabilities of the
active states:

    ū_j = (the sum of ω_s over the active states) / Z_j

A gene with no listed promoter states is not regulated: its control variable is `1.0`.

### Worked example
Gene Y (`j = 2`) at the state `x = [0.0, 0.0, 0.0, 1200.0, 0.0, 600.0]` (protein X at
1200 nmol/gDW and protein Z at 600 nmol/gDW) and `I = 10.0` µM:

    f_X = 0.36     half of protein X is bound to the inducer at 10 µM, so a = 600 nmol/gDW, and
                   with K = 800 and n = 2, f_X = 600² / (800² + 600²) = 0.36
    f_Z = 0.5      protein Z is at its K = 600 nmol/gDW
    ω_1 = 0.00501                     the leak: no transcription factor is bound
    ω_2 = 0.1185 × 0.36 = 0.04266
    ω_3 = 10.01 × 0.5 = 5.004
    Z_2 = 1 + 0.00501 + 0.04266 + 5.004 = 6.052
    ū_2 = (0.00501 + 0.04266 + 5.004) / 6.052 = 0.835

### Arguments
- `network::MyGeneNetwork`: the network, from `load_network(path)`.
- `j::Int64`: the gene, from 1 to the number of genes `N`.
- `x::Array{Float64,1}`: the state, of length `2N`: the mRNA concentrations, then the protein
  concentrations (nmol/gDW).
- `I::Float64`: the inducer concentration (µM).

### Returns
The control variable, a `Float64` from 0 to 1. Leave `x` unchanged.

### Errors
Throws an `ArgumentError` if `j` is not from 1 to `N`, or if `x` does not have length `2N`.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
x = [0.0, 0.0, 0.0, 1200.0, 0.0, 600.0];
control(network, 2, x, 10.0)       # about 0.835, as in the worked example
control(network, 1, zeros(6), 0.0) # about 0.909: gene X has one active state, with W about 10
control(network, 2, zeros(6), 0.0) # about 0.005: with no protein, only the leak state has weight
```
"""
function control(network::MyGeneNetwork, j::Int64, x::Array{Float64,1}, I::Float64)::Float64

    # Step 1: check the arguments -
    N = length(network.genes); # the number of genes
    # TODO: if j is not from 1 to N, throw an ArgumentError, for example
    #       throw(ArgumentError("gene $j is not a gene of the network (1 to $N)"));
    # TODO: if x does not have length 2N, throw an ArgumentError.

    # Step 2: the listed promoter states of gene j -
    states = network.genes[j].promoter;
    # TODO: if states is empty (length(states) == 0), gene j is not regulated: return 1.0.

    # Step 3: add up the weights. The empty promoter has the weight 1 and is never active -
    active_weight = 0.0; # the sum of ω_s over the active states
    total_weight = 1.0;  # the partition function Z_j: 1 plus the sum of ω_s over the listed states
    for s in 1:length(states)

        # The weight ω_s of state s: start from W, then multiply by the bound fraction of each
        # transcription factor bound in state s -
        weight = states[s].W;
        for b in states[s].bound
            # TODO: multiply weight by bound_fraction(network, b, x, I).
        end

        # TODO: add weight to total_weight.
        # TODO: if states[s].active is true, add weight to active_weight.
    end

    # Step 4: the control variable is active_weight / total_weight.
    # TODO: replace the next line with: return active_weight / total_weight;
    error("Complete control in src/Compute.jl");
end

# ----------------------------------------------------------------------------
# Part 1b: The mRNA and protein balances
# ----------------------------------------------------------------------------

"""
    balances!(dx::Array{Float64,1}, x::Array{Float64,1}, parameters::NamedTuple, t::Float64) -> Nothing

Compute the rate of change of every mRNA and protein concentration at the state `x` and the time
`t` (min), and store the rates in `dx` (nmol/gDW/min). The ODE solver in Part 2 calls this function
many times per step. The `!` at the end of the name is a Julia convention: the function changes
its argument `dx`.

### What a material balance says
A material balance on one species in the cell says:

    rate of change of its concentration = rate of production − rates of removal

Each gene j has two balances: one for its mRNA and one for its protein. Transcription of gene j
produces its mRNA, and translation of that mRNA produces its protein. Degradation removes both,
and dilution removes both as the cells grow and divide. Every rate below is in nmol/gDW/min.

### 1. The transcription rate
Compute the transcription rate of gene j in two parts. First, the transcription rate at a control
variable of 1, which is the rate if the promoter of gene j were always in an active state:

    r_X,j = Vmax_X × G / (τ_X × K_X + (1 + τ_X) × G + O_X)

- `Vmax_X` is the transcription kinetic limit. The loader computes it from the RNA polymerase
  concentration, the elongation rate of RNA polymerase, and the length of the gene.
- `G` is the gene concentration, and `K_X` is the transcription saturation constant.
- `τ_X` is the transcription time constant: the time RNA polymerase takes to start transcription,
  divided by the time it takes to transcribe the gene.
- `O_X` accounts for the competition of the other genes for RNA polymerase. It has one term for
  each other gene i. A subscript j marks a value of gene j, and a subscript i a value of gene i:

      O_X = the sum over i ≠ j of (K_X,j × τ_X,j) / (K_X,i × τ_X,i) × (1 + τ_X,i) × G_i

The three genes in `data/memory.toml` have the same values: Vmax_X = 98.209 nmol/gDW/min,
G = 0.0059305 nmol/gDW (one copy per cell), τ_X = 2.52, and K_X = 0.71786 nmol/gDW. Each ratio in
O_X is then 1, and each gene has two other genes, so O_X = 2 × (1 + 2.52) × 0.0059305 = 0.041751
nmol/gDW, and

    r_X,j = 98.209 × 0.0059305 / (2.52 × 0.71786 + (1 + 2.52) × 0.0059305 + 0.041751)
          = 0.31119 nmol/gDW/min

Second, multiply by the control variable ū_j from Part 1a, the probability that the promoter of
gene j is in an active state:

    transcription rate of gene j = r_X,j × ū_j,   with ū_j = control(network, j, x, I)

The control variable depends on the inducer concentration I at the time t. Step 2 of the
function body computes `I = inducer(pulses, t)` for you, so that the transcription rate follows
the inducer pulses.

### 2. The translation rate
The mRNA of gene j is translated into protein at the rate

    translation rate of gene j = r_L,j = Vmax_L × m_j / (τ_L × K_L + (1 + τ_L) × m_j + O_L)

where m_j is the mRNA concentration of gene j.

- `Vmax_L` is the translation kinetic limit. The loader computes it from the ribosome
  concentration, the elongation rate of ribosomes, the length of the protein, and the number of
  ribosomes on each mRNA.
- `K_L` is the translation saturation constant, and `τ_L` is the translation time constant.
- `O_L` accounts for the competition of the other mRNAs for ribosomes. It has one term for each
  other gene i, with the mRNA concentration m_i of gene i:

      O_L = the sum over i ≠ j of (K_L,j × τ_L,j) / (K_L,i × τ_L,i) × (1 + τ_L,i) × m_i

Unlike O_X, O_L changes with the state, because the mRNA concentrations change. The translation
rate is `0` with no mRNA and rises as m_j rises. The model has no control of translation, so this
is the whole translation rate. For every gene in `data/memory.toml`, Vmax_L = 3558.3 nmol/gDW/min,
τ_L = 0.75, and K_L = 136.39 nmol/gDW.

### 3. The degradation rates and the dilution rates
Degradation removes mRNA and protein at a rate proportional to the concentration. The degradation
rate constant is θ = ln 2 / half-life (1/min), one for mRNA and one for protein:

    mRNA degradation rate = θ_m × m_j          protein degradation rate = θ_p × p_j

The cells grow, so the mRNA and the protein spread over more cell mass, and their concentrations
fall. Dilution removes mRNA and protein at a rate proportional to the concentration, with the
specific growth rate μ = ln 2 / doubling time (1/min):

    mRNA dilution rate = μ × m_j               protein dilution rate = μ × p_j

The loader computes θ_m, θ_p, and μ for you. For `data/memory.toml`, θ_m = ln 2 / 5 = 0.13863,
θ_p = ln 2 / 35 = 0.019804, and μ = ln 2 / 60 = 0.011552, all in 1/min.

### 4. The two balances of gene j
    dm_j/dt = transcription rate − mRNA degradation rate − mRNA dilution rate
            = r_X,j × ū_j − θ_m × m_j − μ × m_j

    dp_j/dt = translation rate − protein degradation rate − protein dilution rate
            = r_L,j − θ_p × p_j − μ × p_j

Store the mRNA balance of gene j in `dx[j]` and the protein balance in `dx[N + j]`. For the memory
motif (N = 3), `dx = [dm_X/dt, dm_Y/dt, dm_Z/dt, dp_X/dt, dp_Y/dt, dp_Z/dt]`.

### Where each value is in the code
    gene j                     gene = network.genes[j]
    gene i, in a loop over i   other = network.genes[i]
    m_j, p_j                   x[j], x[N + j]
    G, Vmax_X, τ_X, K_X        gene.G, gene.Vmax_X, gene.τ_X, gene.K_X
    Vmax_L, τ_L, K_L           gene.Vmax_L, gene.τ_L, gene.K_L
    G_i, τ_X,i, K_X,i          other.G, other.τ_X, other.K_X
    m_i, τ_L,i, K_L,i          x[i], other.τ_L, other.K_L
    θ_m, θ_p, μ                network.θ_m, network.θ_p, network.μ
    ū_j                        control(network, j, x, I)
    I                          inducer(pulses, t)
    network, pulses            parameters.network, parameters.pulses

To type τ in VS Code or the Julia REPL, type `\tau` and press Tab.

### Worked example
Gene Y (j = 2) at the state `x = [1.8, 0.05, 0.1, 1200.0, 50.0, 600.0]`, at 70 min in a run with
the pulse `(60.0, 120.0, 10.0)`, so that I = 10 µM. The mRNA concentration of gene Y is
m_Y = x[2] = 0.05 nmol/gDW, the protein concentration is p_Y = x[5] = 50 nmol/gDW, and the other
mRNA concentrations are m_X = x[1] = 1.8 and m_Z = x[3] = 0.1 nmol/gDW. The control variable is
ū_Y = 0.83477, as in the worked example for `control`. Each number is rounded to five
significant figures:

    r_X                      = 0.31119, as in section 1
    transcription rate       = 0.31119 × 0.83477                          = 0.25977
    mRNA degradation rate    = 0.13863 × 0.05                             = 0.0069315
    mRNA dilution rate       = 0.011552 × 0.05                            = 0.00057762
    dx[2] = dm_Y/dt          = 0.25977 − 0.0069315 − 0.00057762           = 0.25226

    O_L                      = (1 + 0.75) × 1.8 + (1 + 0.75) × 0.1        = 3.325
    translation rate         = 3558.3 × 0.05 / (0.75 × 136.39 + (1 + 0.75) × 0.05 + 3.325)
                             = 177.91 / 105.71                            = 1.6831
    protein degradation rate = 0.019804 × 50                              = 0.99021
    protein dilution rate    = 0.011552 × 50                              = 0.57762
    dx[5] = dp_Y/dt          = 1.6831 − 0.99021 − 0.57762                 = 0.11526

### Arguments
- `dx::Array{Float64,1}`: the rates of change, of length `2N`; the function overwrites it.
- `x::Array{Float64,1}`: the state, of length `2N`: the mRNA concentrations, then the protein
  concentrations (nmol/gDW).
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
parameters = (network = network, pulses = [(60.0, 120.0, 10.0)]);
x = [1.8, 0.05, 0.1, 1200.0, 50.0, 600.0];
dx = zeros(6);
balances!(dx, x, parameters, 70.0);
dx[2] # about 0.25226, as in the worked example
dx[5] # about 0.11526
```
"""
function balances!(dx::Array{Float64,1}, x::Array{Float64,1}, parameters::NamedTuple, t::Float64)::Nothing

    # Step 1: unpack the parameters and check the arguments -
    network = parameters.network; # the gene network, from load_network
    pulses = parameters.pulses;   # the inducer pulses
    N = length(network.genes);    # the number of genes
    # TODO: if x does not have length 2N, throw an ArgumentError, for example
    #       throw(ArgumentError("the state has length $(length(x)); it must have length $(2 * N)"));
    # TODO: if dx does not have length 2N, throw an ArgumentError.

    # Step 2: the inducer concentration at time t (µM) -
    I = inducer(pulses, t);

    # Step 3: the degradation rate constants and the specific growth rate (1/min) -
    θ_m = network.θ_m; # the mRNA degradation rate constant, ln 2 / mRNA half-life
    θ_p = network.θ_p; # the protein degradation rate constant, ln 2 / protein half-life
    μ = network.μ;     # the specific growth rate, ln 2 / doubling time

    # Step 4: two balances for each gene j, one for its mRNA and one for its protein -
    for j in 1:N

        # Step 4a: the parameters and the concentrations of gene j -
        gene = network.genes[j];
        m = x[j];     # the mRNA concentration of gene j (nmol/gDW)
        p = x[N + j]; # the protein concentration of gene j (nmol/gDW)

        # Step 4b: the competition of the other genes for RNA polymerase (nmol/gDW). Add one term
        # for each other gene i:
        #     O_X = the sum over i ≠ j of (K_X,j × τ_X,j) / (K_X,i × τ_X,i) × (1 + τ_X,i) × G_i
        # with K_X,j = gene.K_X and τ_X,j = gene.τ_X for gene j, and K_X,i = other.K_X,
        # τ_X,i = other.τ_X, and G_i = other.G for gene i.
        O_X = 0.0;
        for i in 1:N
            if i != j
                other = network.genes[i]; # gene i
                # TODO: O_X = O_X + (the term of gene i in the formula above);
            end
        end

        # Step 4c: the transcription rate of gene j at a control variable of 1 (nmol/gDW/min):
        #     r_X = Vmax_X × G / (τ_X × K_X + (1 + τ_X) × G + O_X)
        # with Vmax_X = gene.Vmax_X, G = gene.G, τ_X = gene.τ_X, K_X = gene.K_X, and O_X from 4b.
        # TODO: compute r_X.

        # Step 4d: the control variable of gene j at the state x and the inducer concentration I,
        # from your function in Part 1a:
        #     u = control(network, j, x, I)
        # TODO: compute u.

        # Step 4e: the transcription rate of gene j (nmol/gDW/min): the transcription rate at a
        # control variable of 1, times the control variable:
        #     transcription_rate = r_X × u
        # TODO: compute transcription_rate.

        # Step 4f: the competition of the other mRNAs for ribosomes (nmol/gDW). Add one term for
        # each other gene i:
        #     O_L = the sum over i ≠ j of (K_L,j × τ_L,j) / (K_L,i × τ_L,i) × (1 + τ_L,i) × m_i
        # with K_L,j = gene.K_L and τ_L,j = gene.τ_L for gene j, and K_L,i = other.K_L,
        # τ_L,i = other.τ_L, and m_i = x[i], the mRNA concentration of gene i.
        O_L = 0.0;
        for i in 1:N
            if i != j
                other = network.genes[i]; # gene i
                # TODO: O_L = O_L + (the term of gene i in the formula above);
            end
        end

        # Step 4g: the translation rate of gene j (nmol/gDW/min). It is 0 with no mRNA and rises as
        # the mRNA concentration m rises:
        #     translation_rate = Vmax_L × m / (τ_L × K_L + (1 + τ_L) × m + O_L)
        # with Vmax_L = gene.Vmax_L, τ_L = gene.τ_L, K_L = gene.K_L, and O_L from 4f.
        # TODO: compute translation_rate.

        # Step 4h: the degradation rates and the dilution rates of the mRNA and the protein of
        # gene j (nmol/gDW/min):
        #     mRNA_degradation_rate = θ_m × m        protein_degradation_rate = θ_p × p
        #     mRNA_dilution_rate = μ × m             protein_dilution_rate = μ × p
        # TODO: compute mRNA_degradation_rate, mRNA_dilution_rate, protein_degradation_rate,
        #       and protein_dilution_rate.

        # Step 4i: the mRNA balance of gene j. Transcription produces mRNA; degradation and
        # dilution remove it:
        #     dm_j/dt = transcription_rate − mRNA_degradation_rate − mRNA_dilution_rate
        # TODO: store dm_j/dt in dx[j].

        # Step 4j: the protein balance of gene j. Translation produces protein; degradation and
        # dilution remove it:
        #     dp_j/dt = translation_rate − protein_degradation_rate − protein_dilution_rate
        # TODO: store dp_j/dt in dx[N + j].
    end

    # Step 5: the rates are in dx, so there is nothing to return.
    # TODO: replace the next line with: return nothing;
    error("Complete balances! in src/Compute.jl");
end

# ----------------------------------------------------------------------------
# Part 2: The ODE solver
# ----------------------------------------------------------------------------

"""
    simulate(network::MyGeneNetwork, x0::Array{Float64,1}, pulses::Array{Tuple{Float64,Float64,Float64},1}; T::Float64) -> NamedTuple

Solve the material balances from time 0 to the final time `T` (min), starting from the state
`x0`, while the inducer concentration follows `pulses`. Return the state at every minute.

Each pulse is a tuple `(start, stop, concentration)`: the inducer concentration is
`concentration` (µM) from `start` up to, but not including, `stop` (min). For example,
`[(60.0, 120.0, 5.0)]` is one pulse of 5 µM from 60 to 120 min. With no pulses, pass
`Tuple{Float64,Float64,Float64}[]`.

### What the ODE solver does
`balances!` gives the rate of change of every concentration at one state and one time. The ODE
solver starts from `x0` at time 0 and takes many small steps forward in time. At each step, it
calls `balances!` to get the rates and uses them to move the state forward. The solver picks the
size of each step itself, so that the error of each step stays below the tolerances you give it.

### Calling the solver
You call the solver with two lines:

    problem = ODEProblem(balances!, copy(x0), (0.0, T), parameters);
    solution = solve(problem, Tsit5(); tstops = stops, saveat = 1.0, reltol = 1e-8, abstol = 1e-12);

The first line describes the problem:
- `balances!`: the function that computes the rates.
- `copy(x0)`: the state at time 0. Pass a copy, so that the solver does not change `x0`.
- `(0.0, T)`: the start time and the final time (min).
- `parameters`: the named tuple `(network = network, pulses = pulses)`, which the solver passes to
  `balances!` as its third argument.

The second line solves the problem:
- `Tsit5()`: the method that the solver uses for each step.
- `tstops = stops`: a list of times that the solver must land on exactly. The inducer
  concentration jumps at the start and at the stop of each pulse, and a solver that steps over a
  jump can miss a short pulse. So `stops` lists the start and the stop of each pulse, if that time
  is between 0 and `T`.
- `saveat = 1.0`: save the state every minute.
- `reltol = 1e-8` and `abstol = 1e-12`: the largest error allowed in each step, relative to each
  concentration (`reltol`) and in nmol/gDW (`abstol`).

### What the solver returns
- `solution.t`: the saved times, `0.0, 1.0, 2.0, ..., T`.
- `solution.u`: the saved states. `solution.u[k]` is the state at the time `solution.t[k]`, a
  vector of length `2N` in the same order as `x0`.
- `solution.retcode`: equal to `ReturnCode.Success` if the solver reached `T`.

`simulate` copies these into two arrays: `t`, the times, and `x`, with one column per time.
Column k of `x` is the state at the time `t[k]`, so `x[5, :]` is protein Y at every minute.

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
length `2N`, or if a pulse stops at or before it starts or has a negative concentration. Throws an
`ErrorException` if the solver stops before `T`.

### Example
```julia
network = load_network(joinpath("data", "memory.toml"));
run = simulate(network, zeros(6), [(60.0, 120.0, 5.0)]; T = 900.0);
run.t[end]    # 900.0
run.x[5, end] # protein Y at 900 min (nmol/gDW)
```
"""
function simulate(network::MyGeneNetwork, x0::Array{Float64,1}, pulses::Array{Tuple{Float64,Float64,Float64},1};
    T::Float64)::NamedTuple

    # Step 1: check the arguments -
    N = length(network.genes); # the number of genes

    # T must be a positive whole number of minutes, such as 900.0 -
    # TODO: throw an ArgumentError if T <= 0.0, if T != round(T), or if isfinite(T) == false.

    # x0 must have 2N entries: the N mRNA concentrations, then the N protein concentrations -
    # TODO: throw an ArgumentError if length(x0) != 2 * N.

    # The pulses: inducer throws an ArgumentError if a pulse is not valid -
    inducer(pulses, 0.0);

    # Step 2: the list of times that the solver must land on exactly: the start and the stop of
    # each pulse, if that time is between 0 and T -
    stops = Float64[]; # an empty list of times
    for k in 1:length(pulses)
        (start, stop, concentration) = pulses[k]; # pulse k
        # TODO: if start > 0.0 and start < T, add start to the list with push!(stops, start);
        # TODO: if stop > 0.0 and stop < T, add stop to the list with push!(stops, stop);
    end

    # Step 3: describe the problem: the function that computes the rates, the state at time 0, the
    # time span, and the parameters that balances! receives -
    parameters = (network = network, pulses = pulses);
    # TODO: problem = ODEProblem(balances!, copy(x0), (0.0, T), parameters);

    # Step 4: solve the problem with the Tsit5 method and the settings in the description -
    # TODO: solution = solve(problem, Tsit5(); tstops = stops, saveat = 1.0,
    #                        reltol = 1e-8, abstol = 1e-12);

    # Step 5: if the solver did not reach T, stop with an error -
    # TODO: if solution.retcode != ReturnCode.Success, call error("the ODE solver stopped early").

    # Step 6: make the arrays to return: t for the times, and x with one column per time -
    # TODO: number_of_times = length(solution.t);
    # TODO: t = zeros(number_of_times);
    # TODO: x = zeros(2 * N, number_of_times);

    # Step 7: copy the solution into t and x, one saved time at a time. solution.t[k] is the k-th
    # saved time, and solution.u[k] is the state at that time, a vector with 2N entries. The state
    # at the k-th saved time goes into column k of x, so entry i of solution.u[k] goes into x[i, k].
    # For example, x[5, 3] is protein Y at the third saved time, 2 min -
    # TODO: uncomment the loops below (delete the # at the start of each line), then replace each
    #       ... with the value that the comment on its line describes.
    # for k in 1:number_of_times          # each saved time
    #     t[k] = ...;                     # the k-th saved time
    #     for i in 1:(2 * N)              # each entry of the state
    #         x[i, k] = ...;              # entry i of the state at the k-th saved time
    #     end
    # end

    # Step 8: return the times and the states.
    # TODO: replace the next line with: return (t = t, x = x);
    error("Complete simulate in src/Compute.jl");
end
