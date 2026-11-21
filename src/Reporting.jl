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
const MAIN_STEP = 0.1;        # the time step of every report run (min)
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
    (A, B) = assignment.state_matrices(memory);
    assignment.control(memory, 2, zeros(6), 0.0);
    assignment.discretize(A, B, 1.0);
    check_run(assignment.simulate(memory, zeros(6), NO_PULSES; Δt = 1.0, T = 2.0), 6, 3);
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
    @printf("Shortest remembered pulse at %.0f µM and Δt = %g min: %d min.\n", MAIN_LEVEL, MAIN_STEP, shortest[k]);
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
    balances(assignment, network, A, B, x, I) -> Array{Float64,1}

Return the right-hand side `A x + B u` of the mRNA and protein balances at the state `x` and the
inducer concentration `I`, with the control `u` from the student's `control`.
"""
function balances(assignment, network, A, B, x, I)

    # The control of each gene at this state -
    N = length(network.genes);
    u = zeros(N);
    for j in 1:N
        u[j] = assignment.control(network, j, x, I);
    end

    # The rate of change of each state -
    return A * x + B * u;
end


"""
    jacobian(assignment, network, x, I) -> Array{Float64,2}

Return the Jacobian of the balances at the state `x` and the inducer concentration `I`, by central
differences: column `c` is the change in the balances per unit change in state `c`.
"""
function jacobian(assignment, network, x, I)

    # Setup -
    (A, B) = assignment.state_matrices(network);
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
        J[:, c] = (balances(assignment, network, A, B, up, I) - balances(assignment, network, A, B, down, I)) / (2 * step);
    end

    return J;
end


"""
    line_state(network, p) -> Array{Float64,1}

Return the state of the memory motif with no inducer in which proteins Y and Z are both at the
level `p` (nmol/gDW), gene X is fully on, and each mRNA has the level that makes its protein
level steady.
"""
function line_state(network, p)

    # Loss rates -
    γ_m = network.θ_m + network.μ;
    γ_p = network.θ_p + network.μ;

    # Gene X is always fully on; the mRNA of Y and Z sets their protein level p -
    m_X = network.genes[1].Vmax_X / γ_m;
    p_X = network.genes[1].k_L * m_X / γ_p;
    m_Y = γ_p * p / network.genes[2].k_L;
    m_Z = γ_p * p / network.genes[3].k_L;

    return [m_X, m_Y, m_Z, p_X, p, p];
end


"""
    line_residual(assignment, network, p) -> Float64

Return the protein level of Y that its control sustains at the state `line_state(network, p)`,
minus `p`. A zero of this function is a steady state of the memory motif with `p_Y = p_Z = p`.
"""
function line_residual(assignment, network, p)

    # The protein level of Y when gene Y is fully on -
    γ_m = network.θ_m + network.μ;
    γ_p = network.θ_p + network.μ;
    full = network.genes[2].k_L * network.genes[2].Vmax_X / (γ_m * γ_p);

    # The level the control sustains, minus p -
    return full * assignment.control(network, 2, line_state(network, p), 0.0) - p;
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

    # The protein level of Y when gene Y is fully on -
    γ_m = network.θ_m + network.μ;
    γ_p = network.θ_p + network.μ;
    full = network.genes[2].k_L * network.genes[2].Vmax_X / (γ_m * γ_p);

    # The grid: 4000 levels of p_Y, evenly spaced on a log scale, from a millionth of the
    # fully on level to twice that level -
    points = 4000;
    log_levels = range(log(1e-6 * full), log(2.0 * full); length = points);
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
            x = line_state(network, (a + b) / 2);

            # The largest real part of the eigenvalues of the Jacobian at this steady state -
            J = jacobian(assignment, network, x, 0.0);
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
    run_pulse(assignment, network, x0, duration, level; Δt = MAIN_STEP) -> NamedTuple

Simulate one inducer pulse of `level` µM that starts at 60 min and lasts `duration` min, from
`x0` to 900 min, with the student's `simulate` and the time step `Δt` (min).
"""
function run_pulse(assignment, network, x0, duration, level; Δt = MAIN_STEP)

    # The pulse, then the run from a copy of x0 -
    pulses = [(PULSE_START, PULSE_START + duration, level)];
    result = assignment.simulate(network, copy(x0), pulses; Δt = Δt, T = FINAL_TIME);

    # Check that the run has one time and one state per step, plus the start -
    return check_run(result, length(x0), round(Int64, FINAL_TIME / Δt) + 1);
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
pulse at 5 µM in the memory motif, from 0 to 900 min, starting at the low steady state. The runs
take steps of `MAIN_STEP`; the table keeps the state every minute, and the plot uses every step.
"""
function write_memory(directory, assignment, memory, no_loop, low, middle)

    # The three runs -
    long = run_pulse(assignment, memory, low, LONG_PULSE, MAIN_LEVEL);
    cut = run_pulse(assignment, no_loop, low, LONG_PULSE, MAIN_LEVEL);
    short = run_pulse(assignment, memory, low, SHORT_PULSE, MAIN_LEVEL);
    pulses = [(PULSE_START, PULSE_START + LONG_PULSE, MAIN_LEVEL)];

    # The table, with one row per minute -
    per_minute = round(Int64, 1.0 / MAIN_STEP); # steps in one minute
    io = open(joinpath(directory, "memory.csv"), "w");
    println(io, "time_min,long_pulse_inducer_uM,pX_nmol_per_gDW,pY_nmol_per_gDW,pZ_nmol_per_gDW,pY_no_loop_nmol_per_gDW,pZ_no_loop_nmol_per_gDW,pY_short_pulse_nmol_per_gDW,pZ_short_pulse_nmol_per_gDW");
    for k in 1:per_minute:length(long.t)
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
with steps of `MAIN_STEP` and find the shortest one that is remembered: protein Y is above its
middle steady state at 900 min. Return the shortest remembered pulse at each level (`0` if none).
"""
function write_pulse_threshold(directory, assignment, memory, low, middle)

    # The header -
    io = open(joinpath(directory, "pulse-threshold.csv"), "w");
    println(io, "inducer_uM,shortest_remembered_min,pY_at_end_of_shortest_nmol_per_gDW,minutes_after_end_to_cross_middle,longest_forgotten_min,peak_pY_after_longest_forgotten_nmol_per_gDW,middle_steady_state_pY_nmol_per_gDW");

    # One row per inducer level -
    shortest = zeros(Int64, length(LEVELS));
    for i in 1:length(LEVELS)
        level = LEVELS[i];

        # Try longer and longer pulses until one is remembered -
        found = 0;
        remembered = nothing;
        d = 1;
        while found == 0 && d <= LONGEST_SCAN
            result = run_pulse(assignment, memory, low, Float64(d), level);
            if result.x[5, end] > middle
                found = d;
                remembered = result;
            end
            d = d + 1;
        end
        shortest[i] = found;
        if found == 0
            @printf(io, "%.0f,none,,,,,%.5f\n", level, middle);
            continue;
        end

        # Protein Y at the end of the shortest remembered pulse, and when it crosses the middle steady state -
        start = round(Int64, PULSE_START / MAIN_STEP) + 1;           # the column of the time the pulse starts
        stop = round(Int64, (PULSE_START + found) / MAIN_STEP) + 1;  # the column of the time the pulse stops
        crossing = stop;
        while remembered.x[5, crossing] <= middle
            crossing = crossing + 1;
        end
        if crossing == stop
            crossing = start;
            while remembered.x[5, crossing] <= middle # crossed during the pulse: find when
                crossing = crossing + 1;
            end
        end

        # The peak of protein Y after the longest forgotten pulse, one minute shorter -
        peak = 0.0;
        if found > 1
            forgotten = run_pulse(assignment, memory, low, Float64(found - 1), level);
            forgotten_stop = round(Int64, (PULSE_START + found - 1) / MAIN_STEP) + 1;
            peak = maximum(forgotten.x[5, forgotten_stop:end]);
        end
        @printf(io, "%.0f,%d,%.5f,%.1f,%d,%.5f,%.5f\n", level, found, remembered.x[5, stop],
            remembered.t[crossing] - remembered.t[stop], found - 1, peak, middle);
    end
    close(io);

    return shortest;
end

end
