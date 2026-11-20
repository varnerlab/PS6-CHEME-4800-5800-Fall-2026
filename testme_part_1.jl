# Run: julia --startup-file=no --project=. testme_part_1.jl
# The flag --solution is for the instructor: it tests solution/src/Compute.jl instead of yours.
# Each labeled case checks one requirement of control or discretize.


# ----------------------------------------------------------------------------
# Setup | Load your functions and the test helpers
# ----------------------------------------------------------------------------

using Test;          # record one scored assertion per labeled case
using LinearAlgebra; # eigvals and the identity matrix I

# Load the code into this script's module; the checker runs each test script in its own module -
Base.include(@__MODULE__, joinpath(@__DIR__, "Include.jl"));
Base.include(@__MODULE__, joinpath(@__DIR__, "test_support.jl"));

using .PS6Checks; # small networks and stored results

# Each case counts as one test, even when your function throws an error in it.
memory = memory_network();
toy = toy_network(TOY_GENES);
(A, B) = state_matrices(memory);
γ_m = memory.θ_m + memory.μ; # mRNA loss rate (1/min)
γ_p = memory.θ_p + memory.μ; # protein loss rate (1/min)

# ----------------------------------------------------------------------------
# Test cases | One scored assertion per labeled case
# ----------------------------------------------------------------------------

checks = [
    # The control of one gene -
    "control returns a Float64" => () ->
        control(memory, 2, copy(LOW_STATE), 0.0) isa Float64,

    "no promoter states gives 1" => () ->
        control(memory, 1, copy(LOW_STATE), 5.0) == 1.0 && control(toy, 1, ones(8), 0.0) == 1.0,

    "leak only" => () -> begin
        u = control(memory, 2, zeros(6), 50.0); # no protein, so only the leak state has weight

        isapprox(u, 0.005 / 1.005; rtol = 1e-12)
    end,

    "inducer switches X on" => () ->
        control(memory, 2, copy(LOW_STATE), 50.0) > control(memory, 2, copy(LOW_STATE), 0.0),

    "partner raises Y" => () -> begin
        x = copy(LOW_STATE);
        x[6] = 3.0; # protein Z at its half-binding level

        control(memory, 2, x, 0.0) > control(memory, 2, copy(LOW_STATE), 0.0)
    end,

    "control between 0 and 1" => () -> begin
        values = zeros(20);
        for k in 1:20
            x = [0.001 * k, 0.0005 * k, 0.0002 * k, 0.6 * k, 0.05 * k^2, 0.2 * k];
            values[k] = control(memory, 2 + (k % 2), x, 2.5 * (k - 1));
        end

        all_between(values, 0.0, 1.0)
    end,

    "Y and Z are symmetric" => () -> begin
        x = [0.01, 0.002, 0.003, 11.0, 0.7, 2.2];
        swapped = [0.01, 0.003, 0.002, 11.0, 2.2, 0.7]; # Y and Z exchanged

        isapprox(control(memory, 2, x, 5.0), control(memory, 3, swapped, 5.0); rtol = 1e-14)
    end,

    "repressor lowers control" => () -> begin
        x = zeros(8);
        x[5] = 2.0; # protein A, which binds in the inactive state of gene D: f_A = 2/3

        isapprox(control(toy, 4, x, 0.0), 0.5 / (1.0 + 0.5 + 4.0 * (2.0 / 3.0)); rtol = 1e-12)
    end,

    "AND state multiplies" => () -> begin
        x = zeros(8);
        x[5] = 1.0; # protein A at its K: f_A = 1/2
        x[6] = 2.0; # protein B at its K: f_B = 1/2

        isapprox(control(toy, 3, x, 0.0), 0.5 / 1.5; rtol = 1e-12)
    end,

    "stored control value" => () ->
        isapprox(control(memory, 2, [0.03, 0.001, 0.002, 11.0, 0.8, 1.5], 5.0), STORED_CONTROL; rtol = 1e-12),

    "bad gene or state rejected" => () ->
        throws_type(() -> control(memory, 0, zeros(6), 0.0), ArgumentError) &&
        throws_type(() -> control(memory, 4, zeros(6), 0.0), ArgumentError) &&
        throws_type(() -> control(memory, 2, zeros(5), 0.0), ArgumentError),

    "state unchanged by control" => () -> begin
        x = [0.03, 0.001, 0.002, 11.0, 0.8, 1.5];
        control(memory, 2, x, 5.0);

        x == [0.03, 0.001, 0.002, 11.0, 0.8, 1.5]
    end,

    # The bilinear discretization -
    "discretize sizes" => () -> begin
        (Ā, B̄) = discretize(A, B, 1.0);

        Ā isa Array{Float64,2} && B̄ isa Array{Float64,2} && size(Ā) == (6, 6) && size(B̄) == (6, 3)
    end,

    "one-state formula" => () -> begin
        (Ā, B̄) = discretize([-0.5;;], [2.0;;], 0.4);

        isapprox(Ā[1, 1], 0.9 / 1.1; rtol = 1e-12) && isapprox(B̄[1, 1], 0.8 / 1.1; rtol = 1e-12)
    end,

    "defining equations" => () -> begin
        (Ā, B̄) = discretize(A, B, 2.0);
        L = I - A; # I - (Δt/2) A with Δt = 2

        isapprox(L * Ā, I + A; atol = 1e-12) && isapprox(L * B̄, 2.0 * B; atol = 1e-12)
    end,

    "eigenvalues of Ā" => () -> begin
        (Ā, B̄) = discretize(A, B, 10.0);
        λ_m = (1 - 5 * γ_m) / (1 + 5 * γ_m); # (1 - Δt γ/2) / (1 + Δt γ/2) with Δt = 10
        λ_p = (1 - 5 * γ_p) / (1 + 5 * γ_p);

        isapprox(sort(real.(eigvals(Ā))), sort([λ_m, λ_m, λ_m, λ_p, λ_p, λ_p]); rtol = 1e-9)
    end,

    "stable at a large step" => () -> begin
        (Ā, B̄) = discretize(A, B, 1000.0);

        maximum(abs.(eigvals(Ā))) < 1.0
    end,

    "bad step rejected" => () ->
        throws_type(() -> discretize(A, B, 0.0), ArgumentError) &&
        throws_type(() -> discretize(A, B, -1.0), ArgumentError),

    "size mismatch rejected" => () ->
        throws_type(() -> discretize(A[:, 1:5], B, 1.0), DimensionMismatch) &&
        throws_type(() -> discretize(A, B[1:5, :], 1.0), DimensionMismatch),

    "matrices unchanged by discretize" => () -> begin
        A_copy = copy(A);
        B_copy = copy(B);
        discretize(A_copy, B_copy, 5.0);

        A_copy == A && B_copy == B
    end,
];


# ----------------------------------------------------------------------------
# Run the 20 scored checks
# ----------------------------------------------------------------------------

@testset "PS6 Part 1 (20 tests)" begin
    for (label, check) in checks # each entry is label => check; check() is true when the case passes
        @testset "$label" begin
            @test check();
        end
    end
end;
