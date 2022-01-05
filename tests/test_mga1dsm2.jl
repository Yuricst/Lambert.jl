"""
Test with mga1dsm
"""

using joptimise
using LinearAlgebra
using Plots
gr()

push!(LOAD_PATH, "../")
using Lambert


# initial and final condition
x01 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
x02 = [1.5, 0.0, 0.0004, 0.0, sqrt(1/1.5), 0.0]
x03 = [0.0, 2.1, 0.002, -sqrt(1/2.1), 0.0, 0.0]
mu =  1.0
μ_bodies = [1.e-4, 1.e-5, 1.e-7]
r_bodies = [4.5e-5, 3e-5, 2e-5]

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

# bounds on problem
tof_max = 4π
t0_bnds = [0.0, 4π]
vinf0_max = 0.05

# get objective function
mga1dsm_problem!, lx, ux, ng, lg, ug = construct_mga1dsm_problem(
	visits,
	mu,
	t0_bnds,
	tof_max,
	vinf0_max,
	r_bodies,
	μ_bodies,
)


## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 200,   # 1500 ~ 2500
    "print_level" => 5,
    "tol" => 1e-6
)
x0 = rand(length(lx))
xopt, fopt, info = minimize(mga1dsm_problem!, x0, ng; 
	lx=lx, 
	ux=ux, 
	lg=lg, 
	ug=ug, 
	solver="ipopt",
	options=ip_options
)

println(xopt)
println(info)

# view solution
visit_nodes, kepler_res, lamb_res, planets = view_mga1dsm_problem(
	xopt, 
	visits,
	mu,
	t0_bnds,
	tof_max,
	vinf0_max,
	r_bodies,
	μ_bodies,
)

# prepare plot
ptraj = plot(
	size=(500,500),
	set_aspect=:equal, 
	legend=false, 
	frame_style=:box, 
)

for res in kepler_res
	plot!(ptraj, res[1,:], res[2,:], c=:green)
end
for lamb in lamb_res
	res = propagate_arc(lamb)
	plot!(ptraj, res[1,:], res[2,:], c=:orange)
end

# inital and final orbit
cs = [:blue, :crimson, :purple]
for (i,planet) in enumerate(planets)
	scatter!(ptraj, [planet[1,1]], [planet[2,1]], marker=:circle, c=cs[i])
	plot!(ptraj, planet[1,:], planet[2,:], linestyle=:dash, c=cs[i])
end
# display plot
display(ptraj)

println("Done!")

