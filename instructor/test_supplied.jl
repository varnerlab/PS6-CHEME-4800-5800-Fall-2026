# Supplied network code checks, separate from the 40 scored student tests.
# Run: julia --startup-file=no --project=. instructor/test_supplied.jl [--solution]
# With the flag, the script also checks the steady states, which need the reference control.
using Test;          # loader and helper assertions
include(joinpath(@__DIR__, "..", "Include.jl"));
include(joinpath(@__DIR__, "..", "src", "Reporting.jl"));

const DATA = joinpath(@__DIR__, "..", "data");
const MEMORY = joinpath(DATA, "memory.toml");
const NO_LOOP = joinpath(DATA, "memory-no-loop.toml");

"""
    write_temp(text) -> String

Write `text` to a new temporary `.toml` file and return its path.
"""
function write_temp(text::String)::String
    path = tempname() * ".toml";
    write(path, text);
    return path;
end

"""
    load_error_message(text) -> String

Load a network from `text` and return the error it throws, as `"Type: message"`, or `""`.
"""
function load_error_message(text::String)::String
    message = "";
    try
        load_network(write_temp(text));
    catch caught
        message = string(nameof(typeof(caught))) * ": " * sprint(showerror, caught);
    end
    return message;
end

# A one-gene network whose promoter has one state bound by the regulator named in REGULATOR -
const ONE_GENE = """
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

  [[genes.promoter]]
  dG = 0.0
  active = true

    [[genes.promoter.bound]]
    regulator = "REGULATOR"
    K = 1.0
    n = 1
""";

# The bound regulator of ONE_GENE with the inducer fields added -
const INDUCER_FIELDS = "n = 1\n    inducer = \"I\"\n    inducer_effect = \"activates\"\n    K_I = 10.0\n    n_I = 1\n";

const RT = 8.314462618e-3 * 310.15; # the thermal energy at 310.15 K (kJ/mol)

@testset "loader" begin

    @testset "loads memory.toml" begin
        network = load_network(MEMORY);
        @test network isa MyGeneNetwork
        @test [g.name for g in network.genes] == ["X", "Y", "Z"]
        @test network.index["Z"] == 3
        @test network.temperature == 310.15
        @test length(network.genes[1].promoter) == 1
        @test length(network.genes[2].promoter) == 3
        @test network.genes[2].promoter[2].bound[1].regulator == "X"
        @test network.genes[2].promoter[2].bound[1].inducer == true
        @test network.genes[2].promoter[3].bound[1].inducer == false
        @test network.genes[2].promoter[3].bound[1].n isa Float64
        for g in network.genes
            @test g.G == 0.005 && g.Vmax_X == 0.04 && g.K_X == 0.005 && g.Vmax_L == 0.05 && g.K_L == 0.05
        end
    end

    @testset "weights from dG" begin
        network = load_network(MEMORY);
        for g in network.genes
            for state in g.promoter
                @test isapprox(state.W, exp(-state.dG / RT); rtol = 1e-15)
            end
        end
        @test network.genes[2].promoter[1].dG == 13.66
        @test isapprox(network.genes[2].promoter[1].W, 0.005006; rtol = 1e-3)
    end

    @testset "loss constants" begin
        network = load_network(MEMORY);
        @test isapprox(network.θ_m, log(2) / 2.5; rtol = 1e-12)
        @test isapprox(network.θ_p, log(2) / 35.0; rtol = 1e-12)
        @test isapprox(network.μ, log(2) / 60.0; rtol = 1e-12)
    end

    @testset "unknown regulator" begin
        message = load_error_message(replace(ONE_GENE, "REGULATOR" => "Q"));
        @test startswith(message, "ArgumentError")
        @test occursin("\"Q\"", message) && occursin("regulator", message)
    end

    @testset "missing field" begin
        message = load_error_message(replace(replace(ONE_GENE, "REGULATOR" => "A"), "K_L = 0.05\n" => ""));
        @test startswith(message, "ArgumentError")
        @test occursin("K_L", message)
    end

    @testset "temperature must be positive" begin
        message = load_error_message(replace(replace(ONE_GENE, "REGULATOR" => "A"), "temperature = 310.15" => "temperature = 0.0"));
        @test startswith(message, "ArgumentError")
        @test occursin("temperature", message)
    end

    @testset "other inducer name" begin
        text = replace(replace(ONE_GENE, "REGULATOR" => "A"), "n = 1\n" => replace(INDUCER_FIELDS, "\"I\"" => "\"J\""));
        message = load_error_message(text);
        @test startswith(message, "ArgumentError")
        @test occursin("inducer", message)
    end

    @testset "inducer needs [inputs.I]" begin
        text = replace(replace(ONE_GENE, "REGULATOR" => "A"), "n = 1\n" => INDUCER_FIELDS);
        message = load_error_message(text);
        @test startswith(message, "ArgumentError")
        @test occursin("inputs.I", message)
        text = replace(text, "[[genes]]" => "[inputs.I]\ninitial = 1.0\n\n[[genes]]");
        message = load_error_message(text);
        @test startswith(message, "ArgumentError")
        @test occursin("initial", message)
        text = replace(text, "initial = 1.0" => "initial = 0.0");
        @test load_error_message(text) == ""
    end

    @testset "gene without promoter states" begin
        text = replace(ONE_GENE, r"  \[\[genes.promoter\]\].*"s => "");
        network = load_network(write_temp(text));
        @test isempty(network.genes[1].promoter)
    end

    @testset "CRLF file loads" begin
        text = replace(read(MEMORY, String), "\n" => "\r\n");
        a = load_network(MEMORY);
        b = load_network(write_temp(text));
        @test [g.name for g in b.genes] == [g.name for g in a.genes]
        @test b.genes[2].promoter[3].W == a.genes[2].promoter[3].W
        @test b.genes[3].promoter[2].bound[1].K == a.genes[3].promoter[2].bound[1].K
        @test b.θ_p == a.θ_p
    end

    @testset "no-loop file" begin
        a = load_network(MEMORY);
        b = load_network(NO_LOOP);
        @test length(b.genes[1].promoter) == 1 && b.genes[1].promoter[1].dG == a.genes[1].promoter[1].dG
        for j in 2:3
            @test length(b.genes[j].promoter) == 2
            for s in 1:2
                @test b.genes[j].promoter[s].dG == a.genes[j].promoter[s].dG
                @test b.genes[j].promoter[s].active == a.genes[j].promoter[s].active
                @test length(b.genes[j].promoter[s].bound) == length(a.genes[j].promoter[s].bound)
            end
            @test b.genes[j].promoter[2].bound[1].K == a.genes[j].promoter[2].bound[1].K
            @test b.genes[j].promoter[2].bound[1].n == a.genes[j].promoter[2].bound[1].n
        end
    end
end

@testset "helpers" begin

    @testset "bound fraction" begin
        network = load_network(MEMORY);
        z_bound = network.genes[2].promoter[3].bound[1]; # Z bound to the promoter of Y
        x_bound = network.genes[2].promoter[2].bound[1]; # X bound to the promoter of Y
        x = zeros(6);
        x[6] = z_bound.K;
        @test bound_fraction(network, z_bound, x, 0.0) == 0.5
        x[4] = 20.0;
        @test bound_fraction(network, x_bound, x, 0.0) == 0.0
        available = 0.5 * x[4]; # half of X carries the inducer at I = K_I, with n_I = 1
        @test isapprox(bound_fraction(network, x_bound, x, x_bound.K_I), available^2 / (x_bound.K^2 + available^2); rtol = 1e-14)
    end

    @testset "negative level binds nothing" begin
        network = load_network(MEMORY);
        b = MyBoundRegulator("Y", 3.0, 2.5, false, 0.0, 0.0);
        x = [0.0, 0.0, 0.0, 1.0, -1e-3, 1.0];
        f = bound_fraction(network, b, x, 0.0);
        @test f === 0.0
    end

    @testset "pulse is half-open" begin
        pulse = [(60.0, 90.0, 5.0)];
        @test inducer(pulse, 90.0) == 0.0
        @test inducer(pulse, 60.0) == 5.0
        @test inducer(pulse, 59.999) == 0.0
        @test inducer([(0.0, 50.0, 2.0), (40.0, 60.0, 3.0)], 45.0) == 5.0
        @test inducer(Tuple{Float64,Float64,Float64}[], 10.0) == 0.0
    end

    @testset "bad pulse" begin
        @test_throws ArgumentError inducer([(90.0, 60.0, 5.0)], 70.0)
        @test_throws ArgumentError inducer([(0.0, 10.0, -1.0)], 5.0)
    end
end

# The steady-state checks call control, so they run only with the reference solution -
if "--solution" in ARGS
@testset "analyses" begin

    @testset "steady states" begin
        network = load_network(MEMORY);
        parameters = (network = network, pulses = Tuple{Float64,Float64,Float64}[]);
        states = MotifReports.steady_states(MemoryMotif, network);
        @test length(states) == 3
        @test [s.largest_real_eigenvalue < 0.0 for s in states] == [true, false, true]
        for s in states
            dx = zeros(6);
            balances!(dx, s.x, parameters, 0.0);
            @test maximum(abs.(dx)) < 1e-12
            @test s.x[5] == s.x[6]
        end
        @test states[1].x[5] < states[2].x[5] < states[3].x[5]
    end

end
end
