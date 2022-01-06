"""
Test with mga1dsm with 1 fly-by, using SPICE ephemeris
"""

using joptimise
using LinearAlgebra
using SPICE
using Plots
gr()

push!(LOAD_PATH, "../")
using Lambert

# modify paths as necessary!
if Sys.iswindows()
    spice_dir = "C:\\Users\\yshimane3\\Documents\\spice"
elseif Sys.islinux()
    spice_dir = "/mnt/c/Users/yurio/Documents/spice"
elseif Sys.isapple()
    spice_dir = "/Users/yuri/Documents/spice"
end
# get spice kernels
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))

# launch window
et0_str = "2022-01-01T00:00:00.00"
etf_str = "2023-01-01T00:00:00.00"
et0 = utc2et(et0_str)
etf = utc2et(et0_str)

locate_x1, cparams = get_spice_locate_function(3, et0_str)
locate_x2, _       = get_spice_locate_function(4, et0_str)
locate_x3, _       = get_spice_locate_function(5, et0_str)

# dynamics constants
mu =  1.0
μ_bodies = [1.e-4, 1.e-5, 1.e-7]
r_bodies = [4.5e-5, 3e-5, 2e-5]

visits = [
	locate_x1, locate_x2, locate_x3
]

# bounds on problem
tof_max = 10*365*86400/cparams.tstar
tf = (etf - et0)/cparams.tstar
t0_bnds = [0.0, tf]
vinf0_max = 5.0/cparams.vstar

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
    "max_iter" => 1000,   # 1500 ~ 2500
    "print_level" => 5,
    "tol" => 1e-5,
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

# view solution
visit_nodes, kepler_res, lamb_res, planets = view_mga1dsm_problem(
	cparams,
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

