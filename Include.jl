# Load the assignment from its own directory, independent of the working directory.
"""
    MemoryMotif

Load the supplied network code, then the three functions in `src/Compute.jl`. The loader,
`load_network(path)`, reads a gene network from a TOML file in `data/`. The helpers
`bound_fraction` and `inducer` return the bound fraction of one transcription factor and the
inducer concentration at a given time. Your functions and the reference solution use the
same supplied code.

The flag `--solution`, passed after a script's filename, is for the instructor: it
loads the instructor's `solution/src/Compute.jl` instead of `src/Compute.jl`. That
folder is not part of the student download. Without the flag, this file always loads
`src/Compute.jl`.

Include this file once per Julia session to make these functions available. After you
edit a source file, restart Julia, or run the scripts from the terminal as the README
shows.
"""
module MemoryMotif

using LinearAlgebra   # eigvals and norm
using OrdinaryDiffEq  # ODEProblem, solve, and the Tsit5 method
using TOML            # TOML.parsefile reads a network file

# Include the supplied network code first, because your functions use it -
include(joinpath(@__DIR__, "src", "Network.jl"));
# If you add helper files, include them here, before Compute.jl. If one helper file
# uses another, include the one it uses first. For example, after creating src/MyHelpers.jl:
# include(joinpath(@__DIR__, "src", "MyHelpers.jl"));

# Choose which Compute.jl to load: the instructor's with --solution, otherwise yours -
if "--solution" in ARGS
    const COMPUTE_PATH = joinpath(@__DIR__, "solution", "src", "Compute.jl");
else
    const COMPUTE_PATH = joinpath(@__DIR__, "src", "Compute.jl");
end
if !isfile(COMPUTE_PATH)
    error("Cannot find $COMPUTE_PATH. The --solution flag needs the instructor's solution folder.");
end
include(COMPUTE_PATH);

export MyBoundRegulator, MyPromoterState, MyGene, MyGeneNetwork;
export load_network, bound_fraction, inducer;
export control, balances!, simulate;

end

using .MemoryMotif; # make the assignment functions available to the calling script
