# ----------------------------------------------------------------------------
# Supplied reports | Run the analyses and write the plots and tables.
# ----------------------------------------------------------------------------

"""Supplied analyses, plots, and tables for the PS6 memory motif."""
module MotifReports

using LinearAlgebra # eigvals
using Plots         # the protein levels over time
using Printf        # tables with a fixed number of digits

# Settings used for every report -
const PULSE_START = 60.0;     # every inducer pulse starts at 60 min
const FINAL_TIME = 900.0;     # every report run ends at 900 min
const MAIN_LEVEL = 5.0;       # the main inducer level (µM)
const LONG_PULSE = 60.0;      # the long pulse in memory.csv lasts 60 min
const SHORT_PULSE = 10.0;     # the short pulse in memory.csv lasts 10 min
const LEVELS = [2.0, 5.0, 10.0, 20.0, 50.0]; # inducer levels in pulse-threshold.csv (µM)
const LONGEST_SCAN = 600;     # the longest pulse tried in pulse-threshold.csv (min)
const NO_PULSES = Tuple{Float64,Float64,Float64}[];

# Every file write_outputs creates; the checker clears exactly these before a run -
const REPORT_FILES = ["steady-states.csv", "memory.csv", "memory.png", "pulse-threshold.csv"];

"""
    write_outputs(assignment, directory) -> NamedTuple

Run the analyses with the student's functions and write every report.

### Arguments
- `assignment`: the loaded `MemoryMotif` module, including the student's three completed
  functions.
- `directory::String`: the output folder; it is created if needed.

### Returns
The named tuple `(middle_pY, shortest_at_main_level)`: the protein level of Y at the middle
steady state, and the shortest remembered pulse (min) at the main inducer level. The function also
writes or overwrites the four files in `REPORT_FILES` and prints a summary.

### Errors
Throws an `ErrorException` if `simulate` returns results of the wrong size, or if the network
does not have three steady states. Errors from the student functions, file writing, or
plotting are passed on.
"""
function write_outputs(assignment, directory)

    # Load the two networks -
    data = joinpath(@__DIR__, "..", "data");
    memory = assignment.load_network(joinpath(data, "memory.toml"));
    no_loop = assignment.load_network(joinpath(data, "memory-no-loop.toml"));

    # Call each student function once on a tiny input, so an unfinished one stops here -
    assignment.control(memory, 2, zeros(6), 0.0);
    assignment.balances!(zeros(6), zeros(6), (network = memory, pulses = NO_PULSES), 0.0);
    check_run(assignment.simulate(memory, zeros(6), NO_PULSES; T = 2.0), 6, 3);
    mkpath(directory);

    # The steady states with no inducer -
    println("Steady states with no inducer.");
    states = steady_states(assignment, memory);
    if length(states) != 3
        error("Found $(length(states)) steady states with no inducer; the memory motif has 3.");
    end
    low = states[1].x;
    middle = states[2].x[5];
    write_steady_states(directory, states);

    # The long pulse in both networks, and the short pulse in the memory motif -
    println("Long and short pulses at ", short_text([MAIN_LEVEL]), " µM.");
    write_memory(directory, assignment, memory, no_loop, low, middle);

    # The shortest remembered pulse at each inducer level -
    println("Shortest remembered pulse at ", short_text(LEVELS), " µM.");
    shortest = write_pulse_threshold(directory, assignment, memory, low, middle);

    # The position of the main inducer level in LEVELS -
    k = 0;
    for i in 1:length(LEVELS)
        if LEVELS[i] == MAIN_LEVEL
            k = i;
        end
    end

    # Print the summary -
    println();
    @printf("Steady states (p_Y, nmol/gDW): low %.4f, middle %.4f, high %.4f.\n", states[1].x[5], middle, states[3].x[5]);
    @printf("Shortest remembered pulse at %.0f µM: %d min.\n", MAIN_LEVEL, shortest[k]);
    println("Reports written to ", directory);

    return (middle_pY = middle, shortest_at_main_level = shortest[k]);
end


"""
    short_text(values) -> String

Return the numbers in `values` as text separated by commas, without trailing zeros: `[0.5, 1.0]`
gives `"0.5, 1"`.
"""
function short_text(values)

    # Write each number with %g, which drops trailing zeros -
    parts = String[];
    for v in values
        push!(parts, @sprintf("%g", v));
    end

    # Join the numbers with commas -
    return join(parts, ", ");
end


"""
    check_run(result, n, points)

Throw an `ErrorException` unless `result`, from the student's `simulate`, has a field `t` with
`points` times and a field `x` with `n` rows and `points` columns. Return `result` when the
sizes are right.
"""
function check_run(result, n, points)

    # Check the times and the states -
    if length(result.t) != points || size(result.x) != (n, points)
        error("simulate returned $(length(result.t)) times and states of size $(size(result.x)); expected $points times and $n × $points states.");
    end

    return result;
end


"""
    rates(assignment, network, x) -> Array{Float64,1}

Return the rate of change of every state at the state `x` with no inducer, from the student's
`balances!` at time 0.
"""
function rates(assignment, network, x)

    # The balances with no inducer pulses -
    dx = zeros(length(x));
    assignment.balances!(dx, copy(x), (network = network, pulses = NO_PULSES), 0.0);

    return dx;
end


"""
    jacobian(assignment, network, x) -> Array{Float64,2}

Return the Jacobian of the balances at the state `x` with no inducer, by central differences:
column `c` is the change in the balances per unit change in state `c`.
"""
function jacobian(assignment, network, x)

    # Setup -
    n = length(x);
    J = zeros(n, n);

    # One central difference per state. The step is a millionth of the state, or a millionth of
    # 1e-4 nmol/gDW for a state near zero -
    for c in 1:n
        step = 1e-6 * max(abs(x[c]), 1e-4);
        up = copy(x);
        down = copy(x);
        up[c] = up[c] + step;
        down[c] = down[c] - step;
        J[:, c] = (rates(assignment, network, up) - rates(assignment, network, down)) / (2 * step);
    end

    return J;
end


"""
    largest_protein_level(network) -> Float64

Return the protein level of gene Y when its mRNA is far above `K_L`, so that translation runs at
`Vmax_L`: no steady state can reach this level (nmol/gDW).
"""
function largest_protein_level(network)

    # Translation at its kinetic limit, balanced by degradation and dilution -
    γ_p = network.θ_p + network.μ;
    return network.genes[2].Vmax_L / γ_p;
end


"""
    line_state(assignment, network, p) -> Array{Float64,1}

Return the state of the memory motif with no inducer in which proteins Y and Z are both at the
level `p` (nmol/gDW), and every mRNA and protein X are steady except the mRNA of Y and Z, which
have the level that makes their protein level steady. `p` must be below
`largest_protein_level(network)`.
"""
function line_state(assignment, network, p)

    # Loss rates -
    γ_m = network.θ_m + network.μ;
    γ_p = network.θ_p + network.μ;

    # Gene X: its promoter state binds no regulator, so its control does not depend on the state -
    X = network.genes[1];
    r_X = X.Vmax_X * X.G / (X.K_X + X.G);
    u_X = assignment.control(network, 1, zeros(6), 0.0);
    m_X = r_X * u_X / γ_m;
    p_X = X.Vmax_L * m_X / (X.K_L + m_X) / γ_p;

    # Genes Y and Z: the mRNA level whose translation balances the loss of protein at level p -
    Y = network.genes[2];
    m = γ_p * p * Y.K_L / (Y.Vmax_L - γ_p * p);

    return [m_X, m, m, p_X, p, p];
end


"""
    line_residual(assignment, network, p) -> Float64

Return the mRNA level of Y that its transcription sustains at the state
`line_state(assignment, network, p)`, minus the mRNA level of Y in that state. A zero of this
function is a steady state of the memory motif with `p_Y = p_Z = p`.
"""
function line_residual(assignment, network, p)

    # The loss rate of mRNA and the transcription rate of gene Y at a control of 1 -
    γ_m = network.θ_m + network.μ;
    Y = network.genes[2];
    r_X = Y.Vmax_X * Y.G / (Y.K_X + Y.G);

    # The level the control sustains, minus the level in the state -
    x = line_state(assignment, network, p);
    return r_X * assignment.control(network, 2, x, 0.0) / γ_m - x[2];
end


"""
    steady_states(assignment, network) -> Array{NamedTuple,1}

Find the steady states of the memory motif with no inducer. Proteins Y and Z play the same role,
so each steady state has `p_Y = p_Z`. The function scans `p_Y` on a logarithmic grid, finds each
sign change of `line_residual`, and narrows it by bisection.

### Returns
One named tuple `(x = x, largest_real_eigenvalue = λ)` per steady state, from the lowest `p_Y`
to the highest: the state and the largest real part of the eigenvalues of the Jacobian (1/min).
A negative value means the steady state is stable.
"""
function steady_states(assignment, network)

    # The grid: 4000 levels of p_Y, evenly spaced on a log scale, from a millionth of the
    # largest protein level to just below it -
    largest = largest_protein_level(network);
    points = 4000;
    log_levels = range(log(1e-6 * largest), log(0.999 * largest); length = points);
    grid = zeros(points);
    for k in 1:points
        grid[k] = exp(log_levels[k]);
    end

    # Bracket each sign change, then bisect -
    states = NamedTuple[];
    for k in 1:(length(grid) - 1)
        a = grid[k];
        b = grid[k + 1];
        f_a = line_residual(assignment, network, a);
        f_b = line_residual(assignment, network, b);
        if sign(f_a) != sign(f_b)

            # Halve the bracket 200 times, more than enough to reach the precision of a Float64 -
            for iteration in 1:200
                c = (a + b) / 2;
                f_c = line_residual(assignment, network, c);
                if sign(f_c) == sign(f_a)
                    a = c;
                    f_a = f_c;
                else
                    b = c;
                end
            end
            x = line_state(assignment, network, (a + b) / 2);

            # The largest real part of the eigenvalues of the Jacobian at this steady state -
            J = jacobian(assignment, network, x);
            eigenvalues = eigvals(J);
            λ = -Inf;
            for i in 1:length(eigenvalues)
                if real(eigenvalues[i]) > λ
                    λ = real(eigenvalues[i]);
                end
            end
            push!(states, (x = x, largest_real_eigenvalue = λ));
        end
    end

    return states;
end


"""
    run_pulse(assignment, network, x0, duration, level) -> NamedTuple

Simulate one inducer pulse of `level` µM that starts at 60 min and lasts `duration` min, from
`x0` to 900 min, with the student's `simulate`.
"""
function run_pulse(assignment, network, x0, duration, level)

    # The pulse, then the run from a copy of x0 -
    pulses = [(PULSE_START, PULSE_START + duration, level)];
    result = assignment.simulate(network, copy(x0), pulses; T = FINAL_TIME);

    # Check that the run has one time and one state per minute, from 0 to 900 min -
    return check_run(result, length(x0), round(Int64, FINAL_TIME) + 1);
end


"""
    write_steady_states(directory, states)

Write `steady-states.csv`: the protein levels and the largest real eigenvalue of each steady state.
"""
function write_steady_states(directory, states)

    # The header -
    names = ["low", "middle", "high"];
    io = open(joinpath(directory, "steady-states.csv"), "w");
    println(io, "state,pX_nmol_per_gDW,pY_nmol_per_gDW,pZ_nmol_per_gDW,largest_real_eigenvalue_per_min");

    # One row per steady state, from low to high -
    for k in 1:3
        x = states[k].x;
        @printf(io, "%s,%.4f,%.5f,%.5f,%.4f\n", names[k], x[4], x[5], x[6], states[k].largest_real_eigenvalue);
    end
    close(io);
end


"""
    write_memory(directory, assignment, memory, no_loop, low, middle)

Write `memory.csv` and `memory.png`: the 60-min pulse at 5 µM in both networks and the 10-min
pulse at 5 µM in the memory motif, from 0 to 900 min, starting at the low steady state. The table
and the plot use the state at every minute.
"""
function write_memory(directory, assignment, memory, no_loop, low, middle)

    # The three runs -
    long = run_pulse(assignment, memory, low, LONG_PULSE, MAIN_LEVEL);
    cut = run_pulse(assignment, no_loop, low, LONG_PULSE, MAIN_LEVEL);
    short = run_pulse(assignment, memory, low, SHORT_PULSE, MAIN_LEVEL);
    pulses = [(PULSE_START, PULSE_START + LONG_PULSE, MAIN_LEVEL)];

    # The table, with one row per minute -
    io = open(joinpath(directory, "memory.csv"), "w");
    println(io, "time_min,long_pulse_inducer_uM,pX_nmol_per_gDW,pY_nmol_per_gDW,pZ_nmol_per_gDW,pY_no_loop_nmol_per_gDW,pZ_no_loop_nmol_per_gDW,pY_short_pulse_nmol_per_gDW,pZ_short_pulse_nmol_per_gDW");
    for k in 1:length(long.t)
        @printf(io, "%.0f,%.1f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", long.t[k], assignment.inducer(pulses, long.t[k]),
            long.x[4, k], long.x[5, k], long.x[6, k], cut.x[5, k], cut.x[6, k], short.x[5, k], short.x[6, k]);
    end
    close(io);

    # The plot of protein Y, on a log scale so that the low state shows; the legend sits below
    # the plot, in two columns, so that it covers no line. Start from the default plot style -
    theme(:default);
    figure = plot(; xlabel = "Time (min)", ylabel = "Protein Y (nmol/gDW)", yscale = :log10,
        title = "Protein Y after a pulse of " * short_text([MAIN_LEVEL]) * " µM inducer",
        legend = :outerbottom, legend_column = 2, size = (900, 560),
        framestyle = :box, titlefontsize = 11, left_margin = 6 * Plots.mm, right_margin = 6 * Plots.mm);
    vspan!(figure, [PULSE_START, PULSE_START + LONG_PULSE]; color = :gray, alpha = 0.15, label = "60-min pulse");
    plot!(figure, long.t, long.x[5, :]; lw = 2, color = 1, label = "Memory motif, 60-min pulse");
    plot!(figure, cut.t, cut.x[5, :]; lw = 2, color = 2, label = "Y–Z loop cut, 60-min pulse");
    plot!(figure, short.t, short.x[5, :]; lw = 2, color = 3, label = "Memory motif, 10-min pulse");
    hline!(figure, [middle]; ls = :dash, lw = 1, color = :black, label = "Middle steady state (memory motif)");
    savefig(figure, joinpath(directory, "memory.png"));
end


"""
    write_pulse_threshold(directory, assignment, memory, low, middle) -> Array{Int64,1}

Write `pulse-threshold.csv`. For each inducer level in `LEVELS`, try pulses of 1, 2, ..., 600 min
and find the shortest one that is remembered: protein Y is above its middle steady state at
900 min. Return the shortest remembered pulse at each level (`0` if none).
"""
function write_pulse_threshold(directory, assignment, memory, low, middle)

    # The header -
    io = open(joinpath(directory, "pulse-threshold.csv"), "w");
    println(io, "inducer_uM,shortest_remembered_min");

    # One row per inducer level -
    shortest = zeros(Int64, length(LEVELS));
    for i in 1:length(LEVELS)
        level = LEVELS[i];

        # Try longer and longer pulses until one is remembered -
        found = 0;
        d = 1;
        while found == 0 && d <= LONGEST_SCAN
            result = run_pulse(assignment, memory, low, Float64(d), level);
            if result.x[5, end] > middle
                found = d;
            end
            d = d + 1;
        end
        shortest[i] = found;

        # Write the row; "none" means no pulse up to LONGEST_SCAN min is remembered -
        if found == 0
            @printf(io, "%.0f,none\n", level);
        else
            @printf(io, "%.0f,%d\n", level, found);
        end
    end
    close(io);

    return shortest;
end

end
