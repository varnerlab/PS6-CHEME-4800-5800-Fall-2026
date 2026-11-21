# Run: julia --startup-file=no --project=. testme_part_2.jl
# The flag --solution is for the instructor: it tests solution/src/Compute.jl instead of yours.
# Each labeled case checks one requirement of simulate.
# simulate calls your control and discretize, so finish Part 1 first. The four cases that check
# the arguments pass on their own, because simulate checks its arguments before any other work.


# ----------------------------------------------------------------------------
# Setup | Load your functions and the test helpers
# ----------------------------------------------------------------------------

using Test;          # record one scored assertion per labeled case
using LinearAlgebra; # the backslash solve for the toy steady state

# Load the code into this script's module; the checker runs each test script in its own module -
Base.include(@__MODULE__, joinpath(@__DIR__, "Include.jl"));
Base.include(@__MODULE__, joinpath(@__DIR__, "test_support.jl"));

using .PS6Checks; # small networks and stored results

# Each case counts as one test, even when your function throws an error in it.
memory = memory_network();
no_loop = no_loop_network();
two = toy_network(TWO_GENES);
long_pulse = [(60.0, 120.0, 5.0)]; # 5 µM from 60 to 120 min

# ----------------------------------------------------------------------------
# Test cases | One scored assertion per labeled case
# ----------------------------------------------------------------------------

checks = [
    # What simulate returns -
    "returns t and x" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 1.0, T = 10.0);

        result.t isa Array{Float64,1} && result.x isa Array{Float64,2}
    end,

    "sizes" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 1.0, T = 900.0);

        length(result.t) == 901 && size(result.x) == (6, 901)
    end,

    "starts at x0" => () -> begin
        x0 = [0.01, 0.0, 0.002, 5.0, 0.0, 1.0];
        result = simulate(memory, x0, NO_PULSES; Δt = 2.0, T = 20.0);

        result.x[:, 1] == x0
    end,

    "times" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 2.5, T = 50.0);

        result.t == 2.5 .* collect(0:20)
    end,

    "T = 1, Δt = 0.1" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 0.1, T = 1.0);

        length(result.t) == 11 && isapprox(result.t[end], 1.0; atol = 1e-12)
    end,

    # The recursion -
    "constant input closed form" => () -> begin
        # With a control of 1, each mRNA moves toward Vmax_X / γ_m by the factor λ each step -
        Δt = 2.0;
        γ_m = two.θ_m + two.μ;
        λ = (1 - Δt * γ_m / 2) / (1 + Δt * γ_m / 2);
        x0 = [0.0, 0.05, 1.0, 0.0];
        result = simulate(two, x0, NO_PULSES; Δt = Δt, T = 100.0);
        ok = true;
        for j in 1:2
            m_star = two.genes[j].Vmax_X / γ_m;
            for k in 0:50
                expected = m_star + λ^k * (x0[j] - m_star);
                ok = ok && isapprox(result.x[j, k + 1], expected; rtol = 1e-10, atol = 1e-15);
            end
        end

        ok
    end,

    "steady state stays" => () -> begin
        (A, B) = state_matrices(two);
        x_star = -(A \ (B * ones(2))); # the steady state with a control of 1
        result = simulate(two, x_star, NO_PULSES; Δt = 5.0, T = 500.0);
        ok = true;
        for k in 1:size(result.x, 2)
            ok = ok && isapprox(result.x[:, k], x_star; rtol = 1e-10);
        end

        ok
    end,

    "low state stays" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 1.0, T = 900.0);

        maximum(abs.(result.x[5, :] .- LOW_STATE[5])) < 1e-9 * LOW_STATE[5]
    end,

    "long pulse switches Y on" => () -> begin
        result = simulate(memory, copy(LOW_STATE), long_pulse; Δt = 1.0, T = 900.0);

        isapprox(result.x[5, end], LONG_PULSE_PY; rtol = 1e-9)
    end,

    "short pulse forgotten" => () -> begin
        result = simulate(memory, copy(LOW_STATE), [(60.0, 70.0, 5.0)]; Δt = 1.0, T = 900.0);

        isapprox(result.x[5, end], SHORT_PULSE_PY; rtol = 1e-9)
    end,

    "no-loop network forgets" => () -> begin
        result = simulate(no_loop, copy(LOW_STATE), long_pulse; Δt = 1.0, T = 900.0);

        isapprox(result.x[5, end], NO_LOOP_PY; rtol = 1e-9)
    end,

    "stored trajectory" => () -> begin
        result = simulate(memory, copy(LOW_STATE), long_pulse; Δt = 1.0, T = 900.0);
        p_Y = [result.x[5, 71], result.x[5, 121], result.x[5, 201]]; # at 70, 120, and 200 min

        isapprox(p_Y, LONG_PULSE_PY_AT_70_120_200; rtol = 1e-9)
    end,

    "pulse edge off the step grid" => () -> begin
        # With 0.1-min steps, the step from 119.9 min sees both pulses and the step from 120 min sees
        # neither, so the runs match. Times found by adding Δt again and again fall just below 120 min -
        a = simulate(memory, copy(LOW_STATE), [(60.0, 120.0, 5.0)]; Δt = 0.1, T = 200.0);
        b = simulate(memory, copy(LOW_STATE), [(60.0, 119.95, 5.0)]; Δt = 0.1, T = 200.0);

        a.x == b.x
    end,

    "inducer from the step start" => () -> begin
        # A pulse that starts at 1 min misses the first step and changes the second -
        with_pulse = simulate(memory, copy(LOW_STATE), [(1.0, 30.0, 50.0)]; Δt = 1.0, T = 5.0);
        without = simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 1.0, T = 5.0);

        with_pulse.x[:, 2] == without.x[:, 2] && with_pulse.x[2, 3] > without.x[2, 3]
    end,

    # The argument checks -
    "bad step rejected" => () ->
        throws_type(() -> simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 0.0, T = 10.0), ArgumentError) &&
        throws_type(() -> simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = -1.0, T = 10.0), ArgumentError),

    "T shorter than a step rejected" => () ->
        throws_type(() -> simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 1.0, T = 0.5), ArgumentError),

    "T not whole steps rejected" => () ->
        throws_type(() -> simulate(memory, copy(LOW_STATE), NO_PULSES; Δt = 1.0, T = 10.5), ArgumentError),

    "wrong x0 length rejected" => () ->
        throws_type(() -> simulate(memory, zeros(5), NO_PULSES; Δt = 1.0, T = 10.0), ArgumentError),

    "x0 unchanged" => () -> begin
        x0 = copy(LOW_STATE);
        simulate(memory, x0, long_pulse; Δt = 1.0, T = 200.0);

        x0 == LOW_STATE
    end,

    "pulses unchanged" => () -> begin
        pulses = [(60.0, 120.0, 5.0), (300.0, 310.0, 50.0)];
        simulate(memory, copy(LOW_STATE), pulses; Δt = 1.0, T = 400.0);

        pulses == [(60.0, 120.0, 5.0), (300.0, 310.0, 50.0)]
    end,
];


# ----------------------------------------------------------------------------
# Run the 20 scored checks
# ----------------------------------------------------------------------------

@testset "PS6 Part 2 (20 tests)" begin
    for (label, check) in checks # each entry is label => check; check() is true when the case passes
        @testset "$label" begin
            @test check();
        end
    end
end;
