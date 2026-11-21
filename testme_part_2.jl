# Run: julia --startup-file=no --project=. testme_part_2.jl
# The flag --solution is for the instructor: it tests solution/src/Compute.jl instead of yours.
# Each labeled case checks one requirement of simulate.
# simulate solves your balances!, which calls your control, so finish Part 1 first. The three cases
# that check the arguments pass on their own, because simulate checks its arguments before any
# other work.


# ----------------------------------------------------------------------------
# Setup | Load your functions and the test helpers
# ----------------------------------------------------------------------------

using Test; # record one scored assertion per labeled case

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
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; T = 10.0);

        result.t isa Array{Float64,1} && result.x isa Array{Float64,2}
    end,

    "sizes" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; T = 900.0);

        length(result.t) == 901 && size(result.x) == (6, 901)
    end,

    "starts at x0" => () -> begin
        x0 = [0.01, 0.0, 0.002, 0.5, 0.0, 0.1];
        result = simulate(memory, x0, NO_PULSES; T = 20.0);

        result.x[:, 1] == x0
    end,

    "one time per minute" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; T = 50.0);

        result.t == collect(0.0:1.0:50.0)
    end,

    # The solution -
    "low state stays" => () -> begin
        result = simulate(memory, copy(LOW_STATE), NO_PULSES; T = 900.0);

        maximum(abs.(result.x[5, :] .- LOW_STATE[5])) < 1e-6 * LOW_STATE[5]
    end,

    "mRNA closed form" => () -> begin
        # With a control of 1, each mRNA approaches r_X / γ_m as exp(-γ_m t) -
        γ_m = two.θ_m + two.μ;
        x0 = [0.0, 0.05, 0.0, 0.0];
        result = simulate(two, x0, NO_PULSES; T = 100.0);
        ok = true;
        for j in 1:2
            g = two.genes[j];
            m_star = g.Vmax_X * g.G / (g.K_X + g.G) / γ_m;
            for k in 0:100
                expected = m_star + (x0[j] - m_star) * exp(-γ_m * k);
                ok = ok && isapprox(result.x[j, k + 1], expected; rtol = 1e-6, atol = 1e-12);
            end
        end

        ok
    end,

    "steady state stays" => () -> begin
        γ_m = two.θ_m + two.μ;
        γ_p = two.θ_p + two.μ;
        x_star = zeros(4);
        for j in 1:2
            g = two.genes[j];
            x_star[j] = g.Vmax_X * g.G / (g.K_X + g.G) / γ_m;                   # mRNA
            x_star[2 + j] = g.Vmax_L * x_star[j] / (g.K_L + x_star[j]) / γ_p; # protein
        end
        result = simulate(two, copy(x_star), NO_PULSES; T = 500.0);
        ok = true;
        for k in 1:size(result.x, 2)
            ok = ok && isapprox(result.x[:, k], x_star; rtol = 1e-8);
        end

        ok
    end,

    "long pulse switches Y on" => () -> begin
        result = simulate(memory, copy(LOW_STATE), long_pulse; T = 900.0);

        isapprox(result.x[5, end], LONG_PULSE_PY; rtol = 1e-6)
    end,

    "short pulse forgotten" => () -> begin
        result = simulate(memory, copy(LOW_STATE), [(60.0, 70.0, 5.0)]; T = 900.0);

        isapprox(result.x[5, end], SHORT_PULSE_PY; rtol = 1e-6)
    end,

    "no-loop network forgets" => () -> begin
        result = simulate(no_loop, copy(LOW_STATE), long_pulse; T = 900.0);

        isapprox(result.x[5, end], NO_LOOP_PY; rtol = 1e-6)
    end,

    "stored trajectory" => () -> begin
        result = simulate(memory, copy(LOW_STATE), long_pulse; T = 900.0);
        p_Y = [result.x[5, 71], result.x[5, 121], result.x[5, 201]]; # at 70, 120, and 200 min

        isapprox(p_Y, LONG_PULSE_PY_AT_70_120_200; rtol = 1e-6)
    end,

    # The pulse edges -
    "half-minute pulse" => () -> begin
        # The solver must stop at both edges of a 0.5-min pulse, or it can step over the pulse -
        result = simulate(memory, copy(LOW_STATE), [(60.25, 60.75, 50.0)]; T = 100.0);

        isapprox(result.x[5, end], EDGE_PULSE_PY; rtol = 1e-6)
    end,

    "inducer starts at the pulse start" => () -> begin
        # A pulse that starts at 1 min leaves the state at 1 min unchanged and raises mRNA Y after it -
        with_pulse = simulate(memory, copy(LOW_STATE), [(1.0, 30.0, 50.0)]; T = 5.0);
        without = simulate(memory, copy(LOW_STATE), NO_PULSES; T = 5.0);

        isapprox(with_pulse.x[:, 2], without.x[:, 2]; rtol = 1e-9) && with_pulse.x[2, 3] > without.x[2, 3] + 1e-6
    end,

    "pulse from time 0" => () -> begin
        result = simulate(memory, copy(LOW_STATE), [(0.0, 30.0, 50.0)]; T = 100.0);

        result.x[5, 31] > 2.0 * LOW_STATE[5]
    end,

    "pulse past the final time" => () -> begin
        result = simulate(memory, copy(LOW_STATE), [(890.0, 950.0, 5.0)]; T = 900.0);

        length(result.t) == 901 && result.x[5, end] > result.x[5, 891]
    end,

    # The argument checks -
    "bad T rejected" => () ->
        throws_type(() -> simulate(memory, copy(LOW_STATE), NO_PULSES; T = 0.0), ArgumentError) &&
        throws_type(() -> simulate(memory, copy(LOW_STATE), NO_PULSES; T = -1.0), ArgumentError),

    "T not whole minutes rejected" => () ->
        throws_type(() -> simulate(memory, copy(LOW_STATE), NO_PULSES; T = 10.5), ArgumentError),

    "wrong x0 length rejected" => () ->
        throws_type(() -> simulate(memory, zeros(5), NO_PULSES; T = 10.0), ArgumentError),

    "x0 unchanged" => () -> begin
        x0 = copy(LOW_STATE);
        simulate(memory, x0, long_pulse; T = 200.0);

        x0 == LOW_STATE
    end,

    "pulses unchanged" => () -> begin
        pulses = [(60.0, 120.0, 5.0), (300.0, 310.0, 50.0)];
        simulate(memory, copy(LOW_STATE), pulses; T = 400.0);

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
