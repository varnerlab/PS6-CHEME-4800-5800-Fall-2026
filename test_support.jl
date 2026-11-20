# ----------------------------------------------------------------------------
# Supplied test support | Small networks and stored results
# ----------------------------------------------------------------------------

"""Supplied test inputs and stored reference results for the PS6 student functions."""
module PS6Checks

using ..MemoryMotif; # implementation loaded by the current test suite

export throws_type, all_between, toy_network, memory_network, no_loop_network;
export TOY_GENES, TWO_GENES, NO_PULSES;
export LOW_STATE, STORED_CONTROL, LONG_PULSE_PY, LONG_PULSE_PY_AT_70_120_200, SHORT_PULSE_PY, NO_LOOP_PY;

# ----------------------------------------------------------------------------
# Stored results | Made once with the reference solution
# (julia --project=. instructor/make_stored_values.jl --solution). A correct solution
# reproduces them to many digits.
# ----------------------------------------------------------------------------

# The low steady state of data/memory.toml with no inducer, [m_X, m_Y, m_Z, p_X, p_Y, p_Z] (nmol/gDW) -
const LOW_STATE = [0.029344417131681514, 0.00014794489931833744, 0.00014794489931833744, 11.229927745853498, 0.05661760198735617, 0.05661760198735617];

# control(memory, 2, [0.03, 0.001, 0.002, 11.0, 0.8, 1.5], 5.0) -
const STORED_CONTROL = 0.5403073972288875;

# simulate(memory, LOW_STATE, [(60.0, 120.0, 5.0)]; Δt = 1.0, T = 900.0): p_Y at 900 min, then at 70, 120, and 200 min -
const LONG_PULSE_PY = 10.185854027605;
const LONG_PULSE_PY_AT_70_120_200 = [0.1774097208440748, 1.7095332559983787, 9.347252821155342];

# simulate(memory, LOW_STATE, [(60.0, 70.0, 5.0)]; Δt = 1.0, T = 900.0): p_Y at 900 min -
const SHORT_PULSE_PY = 0.05661760199135874;

# simulate(memory-no-loop, LOW_STATE, [(60.0, 120.0, 5.0)]; Δt = 1.0, T = 900.0): p_Y at 900 min -
const NO_LOOP_PY = 0.05587028730739187;

# ----------------------------------------------------------------------------
# Test inputs
# ----------------------------------------------------------------------------

# No inducer pulses -
const NO_PULSES = Tuple{Float64,Float64,Float64}[];

# Four genes. A and B have no promoter states. C has one active state that needs both A and B
# bound (an AND gate): its weight is 2 f_A f_B. D has an active leak state with weight 0.5 and an
# inactive state with A bound (a repressor) with weight 4 f_A. The state has length 8.
const TOY_GENES = """
[parameters]
mRNA_half_life = 2.5
protein_half_life = 35.0
doubling_time = 60.0

[[genes]]
name = "A"
Vmax_X = 0.01
k_L = 10.0

[[genes]]
name = "B"
Vmax_X = 0.01
k_L = 10.0

[[genes]]
name = "C"
Vmax_X = 0.01
k_L = 10.0

  [[genes.promoter]]
  W = 2.0
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
Vmax_X = 0.01
k_L = 10.0

  [[genes.promoter]]
  W = 0.5
  active = true

  [[genes.promoter]]
  W = 4.0
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

[[genes]]
name = "A"
Vmax_X = 0.01
k_L = 10.0

[[genes]]
name = "B"
Vmax_X = 0.02
k_L = 5.0
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
