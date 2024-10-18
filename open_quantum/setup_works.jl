#Pkg.activate("../")
Pkg.activate("quantum", shared=true)

using Revise

# Guard against multiple push!
!( "./MyOpenQuantumBase" in LOAD_PATH) && push!(LOAD_PATH, "./MyOpenQuantumBase/")
!( "./MyOpenQuantumTools" in LOAD_PATH) && push!(LOAD_PATH, "./MyOpenQuantumTools/")
