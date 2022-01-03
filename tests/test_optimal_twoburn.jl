"""
Demonstration of simple two-impulse optimization with Lambert solver
"""

using DifferentialEquations
using joptimise
using LinearAlgebra
using Plots
gr()

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

# get objective function
visits = [locate_x1, locate_x2]
rendezvous!, ng, lg, ug = construct_rdv2imp_problem(visits, mu)

# initial guess
t0_guess = 4.0
tof_guess = 3.0
x0 = [t0_guess, tof_guess]
# bounds on variables
lx = [0.0; 0.0]
ux = [5π; 2π]

## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 500,   # 1500 ~ 2500
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

# view result
prop_traj, prop_r1, prop_r2 = view_rdv2imp_problem(xopt, visits)

# prepare plot
ptraj = plot(
	size=(500,500),
	set_aspect=:equal, 
	legend=false, 
	frame_style=:box, 
)

# intial and final position
r1vev = prop_r1[:,1]
r2vec = prop_r2[:,1]
scatter!(ptraj, [r1vec[1]], [r1vec[2]], c=:blue)
scatter!(ptraj, [r2vec[1]], [r2vec[2]], c=:crimson)
plot!(ptraj, [0.0, r1vec[1]], [0.0, r1vec[2]], color=:black)
plot!(ptraj, [0.0, r2vec[1]], [0.0, r2vec[2]], color=:black)

# transfer trajectory
plot!(ptraj, prop_traj[1,:], prop_traj[2,:], linestyle=:dash)
#plot!(ptraj, sol_fwd, vars=(1,2))

# inital and final orbit
plot!(ptraj, prop_r1[1,:], prop_r1[2,:], linestyle=:dash, c=:blue)
plot!(ptraj, prop_r2[1,:], prop_r2[2,:], linestyle=:dash, c=:crimson)

display(ptraj)
println("Done!")



