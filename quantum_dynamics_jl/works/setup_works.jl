Pkg.activate("../")

using Revise

# Guard against multiple push!
!( "./MyOpenQuantumBase" in LOAD_PATH) && push!(LOAD_PATH, "./MyOpenQuantumBase/")
