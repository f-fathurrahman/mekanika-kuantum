using MyOpenQuantumTools
using OrdinaryDiffEq

# XXX: need to simplify this
function my_solve_schrodinger(A::Annealing, tf::Real; tspan = (0, tf), kwargs...)
    
    u0 = build_u0(A.u0, :v) # why need this?
    # probably for interpolations?

    # parameters?
    p = ODEParams(A.H, float(tf), A.annealing_parameter)
    
    update_func! = function (C, u, p, t)
        update_cache!(C, p.L, p, p(t))
    end
    cache = get_cache(A.H)
    diff_op = DiffEqArrayOperator(cache, update_func = update_func!)
    jac_cache = similar(cache)
    ff = ODEFunction(diff_op, jac= update_func!, jac_prototype = jac_cache)

    prob = ODEProblem{true}(ff, u0, float.(tspan), p)
    alg_keyword_warning(;kwargs...)
    solve(prob; kwargs...)
end


# Define Hamiltonian
# s is dimensionless time s=s/tf where tf is final simulation time
H = DenseHamiltonian(
    [(s)->1-s, (s)->s],
    [σx, σz],
    unit=:ħ
)
# The Hamiltonian is sum of two terms
# H = (1-s)*σx + s*σz
# XXX This Hamiltonian is time-dependent

# Final time of simulation
tf = 20

# Initial condition
u0 = PauliVec[1][2]

annealing = Annealing(H, u0)

#sol = solve_schrodinger(annealing, tf, alg=Tsit5(), abstol=1e-6, reltol=1e-6)

# set arguments
A = annealing
alg = Tsit5()
abstol = 1e-6
reltol = 1e-6
tspan = (0, tf)

u0 = build_u0(A.u0, :v) # not needed?
p = ODEParams(A.H, float(tf), A.annealing_parameter)
update_func! = function (C, u, p, t)
    update_cache!(C, p.L, p, p(t))
end
cache = get_cache(A.H)
diff_op = DiffEqArrayOperator(cache, update_func = update_func!)
jac_cache = similar(cache)
ff = ODEFunction(diff_op, jac= update_func!, jac_prototype = jac_cache)

prob = ODEProblem{true}(ff, u0, float.(tspan), p)
sol = solve(prob)

