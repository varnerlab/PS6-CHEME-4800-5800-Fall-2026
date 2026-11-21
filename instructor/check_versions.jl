# Check that every package PS6 shares with the course environment has the course version.
# Run from the assignment root: julia --startup-file=no instructor/check_versions.jl [course Manifest.toml]
# Exits 1 and prints each mismatch; exits 0 when every shared package matches.
using TOML

"""
    manifest_versions(path::String) -> Dict{String,String}

Return the version of every package in the Manifest.toml at `path`, by name. Standard-library
packages without a version are left out.
"""
function manifest_versions(path::String)::Dict{String,String}
    versions = Dict{String,String}()
    deps = TOML.parsefile(path)["deps"]
    for (name, entries) in deps
        entry = entries[1]
        if haskey(entry, "version")
            versions[name] = entry["version"]
        end
    end
    return versions
end

course_manifest = joinpath(@__DIR__, "..", "..", "CHEME-5800-CourseRepository-Fall-2026", "Manifest.toml")
if length(ARGS) >= 1
    course_manifest = ARGS[1]
end
if !isfile(course_manifest)
    println("Course Manifest.toml not found at $course_manifest; skipping the version check.")
    exit(0)
end

ps6 = manifest_versions(joinpath(@__DIR__, "..", "Manifest.toml"))
course = manifest_versions(course_manifest)
mismatches = 0
for name in sort(collect(keys(ps6)))
    if haskey(course, name) && course[name] != ps6[name]
        println("MISMATCH $name: PS6 $(ps6[name]), course $(course[name])")
        global mismatches = mismatches + 1
    end
end
if !haskey(ps6, "OrdinaryDiffEq")
    println("MISSING OrdinaryDiffEq in PS6's Manifest.toml")
    mismatches = mismatches + 1
end
println(mismatches == 0 ? "All shared packages match the course versions." : "$mismatches problem(s).")
exit(mismatches == 0 ? 0 : 1)
