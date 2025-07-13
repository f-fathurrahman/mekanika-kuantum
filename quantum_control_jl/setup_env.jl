Pkg.activate("quantumcontrol", shared=true)

using Revise

!( "./MyQuantumControl" in LOAD_PATH) && push!(LOAD_PATH, "./MyQuantumControl")
!( "./MyQuantumPropagators" in LOAD_PATH) && push!(LOAD_PATH, "./MyQuantumPropagators")
!( "./MyKrotov" in LOAD_PATH) && push!(LOAD_PATH, "./MyKrotov")
