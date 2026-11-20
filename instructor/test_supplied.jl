# Supplied network code checks, separate from the 40 scored student tests.
# Run: julia --startup-file=no --project=. instructor/test_supplied.jl [--solution]
# With the flag, the script also checks the steady states, which need the reference control.
using Test;          # loader and helper assertions
using LinearAlgebra; # eigvals
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

[[genes]]
name = "A"
Vmax_X = 0.01
k_L = 10.0

  [[genes.promoter]]
  W = 1.0
  active = true

    [[genes.promoter.bound]]
    regulator = "REGULATOR"
    K = 1.0
    n = 1
""";

@testset "loader" begin

    @testset "loads memory.toml" begin
        network = load_network(MEMORY);
        @test network isa MyGeneNetwork
        @test [g.name for g in network.genes] == ["X", "Y", "Z"]
        @test network.index["Z"] == 3
        @test isempty(network.genes[1].promoter)
        @test length(network.genes[2].promoter) == 3
        @test network.genes[2].promoter[2].bound[1].regulator == "X"
        @test network.genes[2].promoter[2].bound[1].inducer == true
        @test network.genes[2].promoter[3].bound[1].inducer == false
        @test network.genes[2].promoter[3].bound[1].n isa Float64
        @test network.genes[1].Vmax_X == 8.475e-3
        @test network.genes[3].k_L == 12.0
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
        message = load_error_message(replace(replace(ONE_GENE, "REGULATOR" => "A"), "k_L = 10.0\n" => ""));
        @test startswith(message, "ArgumentError")
        @test occursin("k_L", message)
    end

    @testset "other inducer name" begin
        text = replace(replace(ONE_GENE, "REGULATOR" => "A"), "n = 1\n" => "n = 1\n    inducer = \"J\"\n    inducer_effect = \"activates\"\n    K_I = 1.0\n    n_I = 1\n");
        message = load_error_message(text);
        @test startswith(message, "ArgumentError")
        @test occursin("inducer", message)
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
        for j in 2:3
            @test b.genes[j].promoter[3].W == 0.0
            @test a.genes[j].promoter[3].W > 0.0
            for s in 1:2
                @test b.genes[j].promoter[s].W == a.genes[j].promoter[s].W
                @test length(b.genes[j].promoter[s].bound) == length(a.genes[j].promoter[s].bound)
            end
            @test b.genes[j].promoter[2].bound[1].K == a.genes[j].promoter[2].bound[1].K
            @test b.genes[j].promoter[3].bound[1].regulator == a.genes[j].promoter[3].bound[1].regulator
        end
        @test b.genes[1].Vmax_X == a.genes[1].Vmax_X
    end
end

@testset "helpers" begin

    @testset "A and B" begin
        network = load_network(MEMORY);
        (A, B) = state_matrices(network);
        γ_m = network.θ_m + network.μ;
        γ_p = network.θ_p + network.μ;
        @test size(A) == (6, 6) && size(B) == (6, 3)
        @test A[4, 1] == 12.0 && A[1, 1] == -γ_m && A[5, 5] == -γ_p
        @test B[2, 2] == 8.475e-3 && B[2, 1] == 0.0 && all(B[4:6, :] .== 0.0)
        @test isapprox(sort(eigvals(A)), sort([-γ_m, -γ_m, -γ_m, -γ_p, -γ_p, -γ_p]); rtol = 1e-12)
    end

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
        (A, B) = state_matrices(network);
        states = MotifReports.steady_states(MemoryMotif, network);
        @test length(states) == 3
        @test [s.largest_real_eigenvalue < 0.0 for s in states] == [true, false, true]
        for s in states
            u = [control(network, j, s.x, 0.0) for j in 1:3];
            @test maximum(abs.(A * s.x + B * u)) < 1e-10
            @test s.x[5] == s.x[6]
        end
        @test states[1].x[5] < states[2].x[5] < states[3].x[5]
    end

end
end
