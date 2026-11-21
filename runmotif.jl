# ----------------------------------------------------------------------------
# Run the analyses and write the reports
# ----------------------------------------------------------------------------
# Run from the assignment root:
#   julia --startup-file=no --project=. runmotif.jl
#   julia --startup-file=no --project=. runmotif.jl --solution
#
# The flag --solution is for the instructor: it loads the reference code and writes
# the reports to solution/outputs/.
#
# After you finish Parts 1 and 2, this script runs the analyses described in the
# README section "Run the analyses": the steady states with no inducer, a long and a
# short inducer pulse, and the shortest remembered pulse at five inducer concentrations. It
# writes the plot and tables to outputs/ and overwrites any earlier files with the
# same names.


# 1. Load your functions, the supplied code, and the report code -
ENV["GKSwstype"] = "100"; # render plots to files without opening a GUI window

include(joinpath(@__DIR__, "Include.jl"));
include(joinpath(@__DIR__, "src", "Reporting.jl"));


# 2. Choose the folder for the reports -
if "--solution" in ARGS
    output_directory = joinpath(@__DIR__, "solution", "outputs");
    println("Running the instructor reference solution.");
else
    output_directory = joinpath(@__DIR__, "outputs");
end


# 3. Run the analyses and write the reports -
try
    MotifReports.write_outputs(MemoryMotif, output_directory);
catch caught
    # An unfinished function gets a short message instead of a long error listing -
    if caught isa ErrorException && startswith(caught.msg, "Complete ")
        println(stderr, "Cannot make the reports yet: ", caught.msg, ".");
        println(stderr, "Finish Parts 1 and 2 and pass testme_part_1.jl and testme_part_2.jl first.");
        if !isinteractive() # stop a script run from the terminal; leave an interactive Julia session running
            exit(1);
        end
    else
        rethrow();
    end
end
