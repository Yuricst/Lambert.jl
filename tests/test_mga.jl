"""
Test for MGA Problem
"""

using Plots
using DifferentialEquations
using joptimise

gr()

include("lambert_fast.jl")
include("twobody.jl")
include("mga_problem.jl")


# initial and final condition
x01 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
x02 = [1.5, 0.0, 0.0004, 0.0, sqrt(1/1.5), 0.0]
x02 = [0.0, 2.1, 0.002, -sqrt(1/2.1), 0.0, 0.0]
mu =  1.0
body_mus = [1.e-4, 1.e-5, 1.e-7]
body_radii = [4.5e-5, 3e-5, 2e-5]
h_safes = [0.0, 0.0, 0.0]

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

function locate_x3(epoch::Float64)
	# get position and velocity vector at epoch
	xf = keplerder_nostm(mu, x03, 0.0, epoch, 1.e-14, 20)
	return xf[1:3], xf[4:6]
end


visits = [
	locate_x1, locate_x2, locate_x3
]

# construct problem
mga_problem!, ng = construct_mga_problem(
	visits,
	body_mus, 
	body_radii,
	h_safes, 
	mu, 
)


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


