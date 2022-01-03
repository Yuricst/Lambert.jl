"""
Demonstration of simple two-impulse optimization with Lambert solver
"""

using DifferentialEquations
using joptimise
using LinearAlgebra
using Plots
gr()

include("twobody.jl") # for plotting

push!(LOAD_PATH, "../")
using Lambert


# initial and final condition
x01 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
x02 = [1.5, 0.0, 0.0004, 0.0, sqrt(1/1.5), 0.0]
mu =  1.0
m = 0

function locate_x1(epoch::Float64)
	# get position and velocity vector at epoch
	xf = keplerder_nostm(mu, x01, 0.0, epoch, 1.e-14, 20)
	return xf[1:3], xf[4:6]
end

function locate_x2(epoch::Float64)
	# get position and velocity vector at epoch
	xf = keplerder_nostm(mu, x02, 0.0, epoch, 1.e-14, 20)
	return xf[1:3], xf[4:6]
end


function rendezvous!(g,x)
	# unpack initial vector
	t0, tof = x
	# get initial and final positions
	r1vec, v1vec = locate_x1(t0)
	r2vec, v2vec = locate_x2(t0+tof)
	# solve Lamert
	res = lambert_fast(r1vec, r2vec, tof, m, mu)
	# compute Delta-V costs
	dv1 = norm(res.v1 - v1vec)
	dv2 = norm(v2vec - res.v2)
	# store constraints (if any)
	g[1] = 0.0
	return dv1+dv2
end


# initial guess
t0_guess = 4.0
tof_guess = 3.0
x0 = [t0_guess, tof_guess]
# bounds on variables
lx = [0.0; 0.0]
ux = [5π; 2π]
# bounds on constriants
lg = [0.0]
ug = [0.0]
# number of constraints
ng = 1


## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 2500,   # 1500 ~ 2500
    "print_level" => 5,
    "tol" => 1e-6
)

xopt, fopt, info = minimize(rendezvous!, x0, ng; 
	lx=lx, 
	ux=ux, 
	lg=lg, 
	ug=ug, 
	solver="ipopt",
	options=ip_options
)

println(xopt)
println(info)

# re-construct optimal transfer
t0, tof = xopt

# get initial and final positions
r1vec, v1vec = locate_x1(t0)
r2vec, v2vec = locate_x2(t0+tof)

# solve Lamert
res = lambert_fast(r1vec, r2vec, tof, m, mu)

# construct problem
trajprob_fwd = TwoBodyProblem(
	0.0, 
	tof, 
	vcat(r1vec, res.v1)[:], 
	mu
)
# propagate transfer arc
sol_fwd = propagate(trajprob_fwd, Tsit5())

# propagate initial and final orbit
trajprob_r1 = TwoBodyProblem(
	0.0, 
	get_period(vcat(r1vec, v1vec)[:], mu), 
	vcat(r1vec, v1vec)[:], 
	mu
)
sol_r1 = propagate(trajprob_r1, Tsit5())

trajprob_r2 = TwoBodyProblem(
	0.0, 
	get_period(vcat(r2vec, v2vec)[:], mu), 
	vcat(r2vec, v2vec)[:], 
	mu
)
sol_r2 = propagate(trajprob_r2, Tsit5())


# prepare plot
ptraj = plot(
	size=(500,500),
	set_aspect=:equal, 
	legend=false, 
	frame_style=:box, 
)
# inital and final orbit
plot!(ptraj, sol_r1, vars=(1,2), linestyle=:dash)
plot!(ptraj, sol_r2, vars=(1,2), linestyle=:dash)

# intial and final position
scatter!(ptraj, [r1vec[1]], [r1vec[2]], color=:green)
scatter!(ptraj, [r2vec[1]], [r2vec[2]], color=:red)
plot!(ptraj, [0.0, r1vec[1]], [0.0, r1vec[2]], color=:black)
plot!(ptraj, [0.0, r2vec[1]], [0.0, r2vec[2]], color=:black)

# transfer trajectory
plot!(ptraj, sol_fwd, vars=(1,2))

# scaing plot
#scatter!(ptraj, [-2.0,2.0], [-2.0,2.0], alpha=0.0)

display(ptraj)
println("Done!")



