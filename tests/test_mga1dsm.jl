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
x02 = [0.9, 0.0, 0.7, 0.0, sqrt(1/sqrt(0.8^2+0.7^2)), 0.0]
mu =  1.0

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

# bounds on problem
tof_max = 4π
t0_bnds = [0.0, 4π]
vinf0_max = 0.05

# get objective function
visits = [locate_x1, locate_x2]
mga1dsm_problem!, lx, ux, ng, lg, ug = construct_mga1dsm_problem(
	visits,
	mu,
	t0_bnds,
	tof_max,
	vinf0_max,
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
cs = [:blue, :crimson]
for (i,planet) in enumerate(planets)
	plot!(ptraj, planet[1,:], planet[2,:], linestyle=:dash, c=cs[i])
end
# display plot
display(ptraj)

println("Done!")

