# ----------------------------------------------------------------------------
# Supplied test support | Small networks and stored results
# ----------------------------------------------------------------------------

"""Supplied test inputs and stored reference results for the PS6 student functions."""
module PS6Checks

using ..MemoryMotif; # implementation loaded by the current test suite

export throws_type, all_between, toy_network, memory_network, no_loop_network;
export TOY_GENES, TWO_GENES, NO_PULSES, X_TEST, transcription_limit, translation_limit;
export LOW_STATE, STORED_CONTROL, STORED_BALANCES;
export LONG_PULSE_PY, LONG_PULSE_PY_AT_70_120_200, SHORT_PULSE_PY, NO_LOOP_PY, EDGE_PULSE_PY;

# ----------------------------------------------------------------------------
# Stored results | Made once with the reference solution
# (julia --project=. instructor/make_stored_values.jl --solution). A correct solution
# reproduces them to many digits.
# ----------------------------------------------------------------------------

# The low steady state of data/memory.toml with no inducer, [m_X, m_Y, m_Z, p_X, p_Y, p_Z] (nmol/gDW) -
const LOW_STATE = [1.8838556568621134, 0.010455210963460002, 0.010455210963460002, 2023.8633214824754, 11.232239535036165, 11.232239535036165];

# control(memory, 2, X_TEST, 5.0) -
const STORED_CONTROL = 0.8347663835474028;

# balances! at the state X_TEST and 70 min, with the pulse (60.0, 120.0, 5.0) (nmol/gDW/min) -
const STORED_BALANCES = [0.05764816769745795, 0.25226051201990946, 0.15208473094723338, -5.69718281293968, -7.7155040541297755, -15.431008108259551];

# simulate(memory, LOW_STATE, [(60.0, 120.0, 5.0)]; T = 900.0): p_Y at 900 min, then at 70, 120, and 200 min -
const LONG_PULSE_PY = 1900.547702642292;
const LONG_PULSE_PY_AT_70_120_200 = [25.203101311095597, 140.983746291264, 1536.0930813517123];

# simulate(memory, LOW_STATE, [(60.0, 70.0, 5.0)]; T = 900.0): p_Y at 900 min -
const SHORT_PULSE_PY = 11.232239532939264;

# simulate(memory-no-loop, LOW_STATE, [(60.0, 120.0, 5.0)]; T = 900.0): p_Y at 900 min -
const NO_LOOP_PY = 11.087579886408271;

# simulate(memory, LOW_STATE, [(60.25, 60.75, 50.0)]; T = 100.0): p_Y at 100 min -
const EDGE_PULSE_PY = 12.3847518255152;

# ----------------------------------------------------------------------------
# Test inputs
# ----------------------------------------------------------------------------

# No inducer pulses -
const NO_PULSES = Tuple{Float64,Float64,Float64}[];

# A state away from every steady state of the memory motif (nmol/gDW) -
const X_TEST = [1.5, 0.05, 0.1, 1800.0, 300.0, 600.0];

# Four genes. A and B have no promoter states. C has one active state that needs both A and B
# bound (an AND gate): its weight is W f_A f_B, with W about 2. D has an active leak state with
# W about 0.5 and an inactive state with A bound (a repressor) with weight W f_A, W about 4.
# The state has length 8.
const TOY_GENES = """
[parameters]
mass_of_single_cell = 2.8e-13
volume_of_single_cell = 6.7e-16
doubling_time = 60.0
temperature = 310.15
copies_of_rnap_per_cell = 4600.0
transcription_elongation_rate = 60.0
transcription_initiation_time_constant = 42.0
transcription_saturation_constant = 0.3
copies_of_ribosome_per_cell = 50000.0
translation_elongation_rate = 16.5
translation_initiation_time_constant = 15.0
translation_saturation_factor = 190.0
average_polysome_number = 4.0
mRNA_half_life = 5.0
protein_half_life = 35.0

[[genes]]
name = "A"
gene_copy_number = 1.0
length_of_gene_in_nt = 1000.0
length_of_protein_in_aa = 330.0

[[genes]]
name = "B"
gene_copy_number = 1.0
length_of_gene_in_nt = 1000.0
length_of_protein_in_aa = 330.0

[[genes]]
name = "C"
gene_copy_number = 1.0
length_of_gene_in_nt = 1000.0
length_of_protein_in_aa = 330.0

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
gene_copy_number = 1.0
length_of_gene_in_nt = 1000.0
length_of_protein_in_aa = 330.0

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

# Two genes with no promoter states, so the control variable of each is always 1. Gene B has two
# copies, a longer gene, and a shorter protein, so the two genes have different kinetic limits and
# time constants -
const TWO_GENES = """
[parameters]
mass_of_single_cell = 2.8e-13
volume_of_single_cell = 6.7e-16
doubling_time = 60.0
temperature = 310.15
copies_of_rnap_per_cell = 4600.0
transcription_elongation_rate = 60.0
transcription_initiation_time_constant = 42.0
transcription_saturation_constant = 0.3
copies_of_ribosome_per_cell = 50000.0
translation_elongation_rate = 16.5
translation_initiation_time_constant = 15.0
translation_saturation_factor = 190.0
average_polysome_number = 4.0
mRNA_half_life = 5.0
protein_half_life = 35.0

[[genes]]
name = "A"
gene_copy_number = 1.0
length_of_gene_in_nt = 1000.0
length_of_protein_in_aa = 330.0

[[genes]]
name = "B"
gene_copy_number = 2.0
length_of_gene_in_nt = 1500.0
length_of_protein_in_aa = 200.0
""";

"""
    transcription_limit(network, j) -> Float64

Return the transcription rate of gene `j` at a control variable of 1 (nmol/gDW/min), with the
competition of the other genes for RNA polymerase. An independent copy of the formula, for the tests.
"""
function transcription_limit(network, j)
    gene = network.genes[j]
    O_X = 0.0
    for i in 1:length(network.genes)
        if i != j
            other = network.genes[i]
            O_X += (gene.K_X * gene.τ_X) / (other.K_X * other.τ_X) * (1.0 + other.τ_X) * other.G
        end
    end
    return gene.Vmax_X * gene.G / (gene.τ_X * gene.K_X + (1.0 + gene.τ_X) * gene.G + O_X)
end

"""
    translation_limit(network, j, x) -> Float64

Return the translation rate of gene `j` at the state `x` (nmol/gDW/min), with the competition of
the other mRNAs for ribosomes. An independent copy of the formula, for the tests.
"""
function translation_limit(network, j, x)
    gene = network.genes[j]
    O_L = 0.0
    for i in 1:length(network.genes)
        if i != j
            other = network.genes[i]
            O_L += (gene.K_L * gene.τ_L) / (other.K_L * other.τ_L) * (1.0 + other.τ_L) * x[i]
        end
    end
    return gene.Vmax_L * x[j] / (gene.τ_L * gene.K_L + (1.0 + gene.τ_L) * x[j] + O_L)
end

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
