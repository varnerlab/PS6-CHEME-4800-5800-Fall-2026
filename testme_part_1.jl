# Run: julia --startup-file=no --project=. testme_part_1.jl
# The flag --solution is for the instructor: it tests solution/src/Compute.jl instead of yours.
# Each labeled case checks one requirement of control or balances!.


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
toy = toy_network(TOY_GENES);
two = toy_network(TWO_GENES);
long_pulse = [(60.0, 120.0, 5.0)];           # 5 µM from 60 to 120 min
x_test = [0.03, 0.001, 0.002, 0.8, 0.3, 0.5]; # a state away from every steady state (nmol/gDW)
γ_m = two.θ_m + two.μ;                        # mRNA loss rate (1/min)

# ----------------------------------------------------------------------------
# Test cases | One scored assertion per labeled case
# ----------------------------------------------------------------------------

checks = [
    # The control of one gene -
    "control returns a Float64" => () ->
        control(memory, 2, copy(LOW_STATE), 0.0) isa Float64,

    "no promoter states gives 1" => () ->
        control(toy, 1, ones(8), 0.0) == 1.0 && control(toy, 2, zeros(8), 5.0) == 1.0,

    "states with no bound regulators" => () -> begin
        W_X = memory.genes[1].promoter[1].W;    # gene X: RNA polymerase bound alone
        W_leak = memory.genes[2].promoter[1].W; # gene Y: the leak state
        u_X = control(memory, 1, copy(LOW_STATE), 5.0);
        u_Y = control(memory, 2, zeros(6), 50.0); # no protein, so only the leak state has weight

        isapprox(u_X, W_X / (1.0 + W_X); rtol = 1e-12) && isapprox(u_Y, W_leak / (1.0 + W_leak); rtol = 1e-12)
    end,

    "inducer switches X on" => () ->
        control(memory, 2, copy(LOW_STATE), 50.0) > control(memory, 2, copy(LOW_STATE), 0.0),

    "partner raises Y" => () -> begin
        x = copy(LOW_STATE);
        x[6] = memory.genes[2].promoter[3].bound[1].K; # protein Z at its half-binding level

        control(memory, 2, x, 0.0) > control(memory, 2, copy(LOW_STATE), 0.0)
    end,

    "control between 0 and 1" => () -> begin
        values = zeros(20);
        for k in 1:20
            x = [0.001 * k, 0.0005 * k, 0.0002 * k, 0.05 * k, 0.003 * k^2, 0.01 * k];
            values[k] = control(memory, 2 + (k % 2), x, 2.5 * (k - 1));
        end

        all_between(values, 0.0, 1.0)
    end,

    "Y and Z are symmetric" => () -> begin
        x = [0.01, 0.002, 0.003, 0.9, 0.07, 0.2];
        swapped = [0.01, 0.003, 0.002, 0.9, 0.2, 0.07]; # Y and Z exchanged

        isapprox(control(memory, 2, x, 5.0), control(memory, 3, swapped, 5.0); rtol = 1e-14)
    end,

    "repressor lowers control" => () -> begin
        x = zeros(8);
        x[5] = 2.0; # protein A, which binds in the inactive state of gene D: f_A = 2/3
        W_leak = toy.genes[4].promoter[1].W;
        W_bound = toy.genes[4].promoter[2].W;

        isapprox(control(toy, 4, x, 0.0), W_leak / (1.0 + W_leak + W_bound * (2.0 / 3.0)); rtol = 1e-12)
    end,

    "AND state multiplies" => () -> begin
        x = zeros(8);
        x[5] = 1.0; # protein A at its K: f_A = 1/2
        x[6] = 2.0; # protein B at its K: f_B = 1/2
        W = toy.genes[3].promoter[1].W;

        isapprox(control(toy, 3, x, 0.0), 0.25 * W / (1.0 + 0.25 * W); rtol = 1e-12)
    end,

    "stored control value" => () ->
        isapprox(control(memory, 2, copy(x_test), 5.0), STORED_CONTROL; rtol = 1e-12),

    "bad gene or state rejected" => () ->
        throws_type(() -> control(memory, 0, zeros(6), 0.0), ArgumentError) &&
        throws_type(() -> control(memory, 4, zeros(6), 0.0), ArgumentError) &&
        throws_type(() -> control(memory, 2, zeros(5), 0.0), ArgumentError),

    "state unchanged by control" => () -> begin
        x = copy(x_test);
        control(memory, 2, x, 5.0);

        x == x_test
    end,

    # The mRNA and protein balances -
    "balances! fills dx and returns nothing" => () -> begin
        dx = fill(NaN, 6); # every entry must be overwritten
        result = balances!(dx, copy(x_test), (network = memory, pulses = long_pulse), 70.0);

        isnothing(result) && all_between(dx, -1.0, 1.0)
    end,

    "stored balances" => () -> begin
        dx = zeros(6);
        balances!(dx, copy(x_test), (network = memory, pulses = long_pulse), 70.0);

        isapprox(dx, STORED_BALANCES; rtol = 1e-12)
    end,

    "zero at the low steady state" => () -> begin
        dx = ones(6);
        balances!(dx, copy(LOW_STATE), (network = memory, pulses = NO_PULSES), 0.0);

        maximum(abs.(dx)) < 1e-12
    end,

    "mRNA balance" => () -> begin
        # Genes A and B have no promoter states, so their control is 1 -
        x = [0.01, 0.03, 0.2, 0.5];
        dx = zeros(4);
        balances!(dx, x, (network = two, pulses = NO_PULSES), 0.0);
        ok = true;
        for j in 1:2
            g = two.genes[j];
            expected = g.Vmax_X * g.G / (g.K_X + g.G) - γ_m * x[j];
            ok = ok && isapprox(dx[j], expected; rtol = 1e-12);
        end

        ok
    end,

    "translation saturates" => () -> begin
        A = two.genes[1];
        B = two.genes[2];
        x = [A.K_L, 100.0 * B.K_L, 0.0, 0.0]; # no protein, so no protein is lost
        dx = zeros(4);
        balances!(dx, x, (network = two, pulses = NO_PULSES), 0.0);

        isapprox(dx[3], A.Vmax_L / 2.0; rtol = 1e-12) && isapprox(dx[4], B.Vmax_L * 100.0 / 101.0; rtol = 1e-12)
    end,

    "inducer read at time t" => () -> begin
        # The pulse is on at 70 min and off from 120 min, when the low state is steady again -
        during = zeros(6);
        after = zeros(6);
        balances!(during, copy(LOW_STATE), (network = memory, pulses = long_pulse), 70.0);
        balances!(after, copy(LOW_STATE), (network = memory, pulses = long_pulse), 120.0);

        during[2] > 1e-6 && abs(after[2]) < 1e-12
    end,

    "wrong lengths rejected" => () ->
        throws_type(() -> balances!(zeros(6), zeros(5), (network = memory, pulses = NO_PULSES), 0.0), ArgumentError) &&
        throws_type(() -> balances!(zeros(5), zeros(6), (network = memory, pulses = NO_PULSES), 0.0), ArgumentError),

    "x unchanged by balances!" => () -> begin
        x = copy(x_test);
        balances!(zeros(6), x, (network = memory, pulses = long_pulse), 70.0);

        x == x_test
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
