# ----------------------------------------------------------------------------
# Supplied test support | Small networks and stored results
# ----------------------------------------------------------------------------

"""Supplied test inputs and stored reference results for the PS6 student functions."""
module PS6Checks

using ..MemoryMotif; # implementation loaded by the current test suite

export throws_type, all_between, toy_network, memory_network, no_loop_network;
export TOY_GENES, TWO_GENES, NO_PULSES;
export LOW_STATE, STORED_CONTROL, STORED_BALANCES;
export LONG_PULSE_PY, LONG_PULSE_PY_AT_70_120_200, SHORT_PULSE_PY, NO_LOOP_PY, EDGE_PULSE_PY;

# ----------------------------------------------------------------------------
# Stored results | Made once with the reference solution
# (julia --project=. instructor/make_stored_values.jl --solution). A correct solution
# reproduces them to many digits.
# ----------------------------------------------------------------------------

# The low steady state of data/memory.toml with no inducer, [m_X, m_Y, m_Z, p_X, p_Y, p_Z] (nmol/gDW) -
const LOW_STATE = [0.06295896471273761, 0.000352553116905135, 0.000352553116905135, 0.8887448708441775, 0.011164603267003845, 0.011164603267003845];

# control(memory, 2, [0.03, 0.001, 0.002, 0.8, 0.3, 0.5], 5.0) -
const STORED_CONTROL = 0.8345969770388688;

# balances! at the state [0.03, 0.001, 0.002, 0.8, 0.3, 0.5] and 70 min, with the pulse (60.0, 120.0, 5.0) (nmol/gDW/min) -
const STORED_BALANCES = [0.009518922277003667, 0.016403128215544068, 0.012328509145562497, -0.006335326534550403, -0.008426605293593655, -0.013755252161017077];

# simulate(memory, LOW_STATE, [(60.0, 120.0, 5.0)]; T = 900.0): p_Y at 900 min, then at 70, 120, and 200 min -
const LONG_PULSE_PY = 0.8823134618527292;
const LONG_PULSE_PY_AT_70_120_200 = [0.026829437228647085, 0.15329262609970745, 0.7889061459824324];

# simulate(memory, LOW_STATE, [(60.0, 70.0, 5.0)]; T = 900.0): p_Y at 900 min -
const SHORT_PULSE_PY = 0.011164603266713183;

# simulate(memory-no-loop, LOW_STATE, [(60.0, 120.0, 5.0)]; T = 900.0): p_Y at 900 min -
const NO_LOOP_PY = 0.010924356108919604;

# simulate(memory, LOW_STATE, [(60.25, 60.75, 50.0)]; T = 100.0): p_Y at 100 min -
const EDGE_PULSE_PY = 0.012160329197283588;

# ----------------------------------------------------------------------------
# Test inputs
# ----------------------------------------------------------------------------

# No inducer pulses -
const NO_PULSES = Tuple{Float64,Float64,Float64}[];

# Four genes. A and B have no promoter states. C has one active state that needs both A and B
# bound (an AND gate): its weight is W f_A f_B, with W about 2. D has an active leak state with
# W about 0.5 and an inactive state with A bound (a repressor) with weight W f_A, W about 4.
# The state has length 8.
const TOY_GENES = """
[parameters]
mRNA_half_life = 2.5
protein_half_life = 35.0
doubling_time = 60.0
temperature = 310.15

[[genes]]
name = "A"
G = 0.005
Vmax_X = 0.04
K_X = 0.005
Vmax_L = 0.05
K_L = 0.05

[[genes]]
name = "B"
G = 0.005
Vmax_X = 0.04
K_X = 0.005
Vmax_L = 0.05
K_L = 0.05

[[genes]]
name = "C"
G = 0.005
Vmax_X = 0.04
K_X = 0.005
Vmax_L = 0.05
K_L = 0.05

  [[genes.promoter]]
  dG = -1.79
  active = true

    [[genes.promoter.bound]]
    regulator = "A"
    K = 1.0
    n = 1

    [[genes.promoter.bound]]
    regulator = "B"
    K = 2.0
    n = 2

[[genes]]
name = "D"
G = 0.005
Vmax_X = 0.04
K_X = 0.005
Vmax_L = 0.05
K_L = 0.05

  [[genes.promoter]]
  dG = 1.79
  active = true

  [[genes.promoter]]
  dG = -3.57
  active = false

    [[genes.promoter.bound]]
    regulator = "A"
    K = 1.0
    n = 1
""";

# Two genes with no promoter states, so the control of each is always 1 -
const TWO_GENES = """
[parameters]
mRNA_half_life = 2.5
protein_half_life = 35.0
doubling_time = 60.0
temperature = 310.15

[[genes]]
name = "A"
G = 0.005
Vmax_X = 0.04
K_X = 0.005
Vmax_L = 0.05
K_L = 0.05

[[genes]]
name = "B"
G = 0.01
Vmax_X = 0.02
K_X = 0.004
Vmax_L = 0.1
K_L = 0.02
""";

"""
    toy_network(text::String) -> MyGeneNetwork

Write `text` to a temporary TOML file and load it with `load_network`.
"""
function toy_network(text::String)::MyGeneNetwork
    path = tempname() * ".toml";
    write(path, text);
    return load_network(path);
end

"""
    memory_network() -> MyGeneNetwork

Load `data/memory.toml`.
"""
function memory_network()::MyGeneNetwork
    return load_network(joinpath(@__DIR__, "data", "memory.toml"));
end

"""
    no_loop_network() -> MyGeneNetwork

Load `data/memory-no-loop.toml`, the motif with the Y–Z loop cut.
"""
function no_loop_network()::MyGeneNetwork
    return load_network(joinpath(@__DIR__, "data", "memory-no-loop.toml"));
end

"""
    throws_type(f, T) -> Bool

Call `f()` and return `true` if it throws an error of type `T`. Return `false` if it
throws another kind of error or finishes without an error.
"""
function throws_type(f, T)

    # Run f; if it throws, report whether the error has the expected type -
    try
        f();
    catch e
        return e isa T;
    end

    # f finished without an error -
    return false;
end

"""
    all_between(values, low, high) -> Bool

Return `true` if every entry of `values` is at least `low` and at most `high`.
"""
function all_between(values, low, high)

    # Stop at the first entry outside the range; NaN is outside every range -
    for k in 1:length(values)
        if !(values[k] >= low && values[k] <= high)
            return false;
        end
    end

    return true;
end

end
