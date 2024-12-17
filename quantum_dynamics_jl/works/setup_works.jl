Pkg.activate("../")

using Revise

# Guard against multiple push!
!( "./MyQuantumDynamics" in LOAD_PATH) && push!(LOAD_PATH, "./MyQuantumDynamics/")
