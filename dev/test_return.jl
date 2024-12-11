"""
Test with mga1dsm with 1 fly-by, using SPICE ephemeris
"""

using joptimise
using LinearAlgebra
using SPICE
using Plots
gr()

#push!(LOAD_PATH, "../")
#using Lambert
include("../src/lambert.jl")

# constants
AU = 1.495978707e8
MU_SUN = 1.3271244004193938e11
MU_EARTH = 4.0350323550225981E+05
MU_VENUS = 3.2485859200000006E+05
MU_SATURN = 3.7940585200000003E+07

# modify paths as necessary!
spice_dir = ENV["SPICE"]

# get spice kernels
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))

# launch window
et0_str = "2030-01-01T00:00:00.00"
etf_str = "2040-01-01T00:00:00.00"
et0 = utc2et(et0_str)
etf = utc2et(et0_str)

locate_x1, cparams = Lambert.get_spice_locate_function(5, et0_str)
locate_x2, _       = Lambert.get_spice_locate_function(5, et0_str)
locate_x3, _       = Lambert.get_spice_locate_function(3, et0_str)

# dynamics constants
mu =  1.0
μ_bodies = [
	MU_SATURN / MU_SUN,
	MU_SATURN / MU_SUN,
    MU_EARTH  / MU_SUN,
]
r_bodies = [4.5e-5, 3e-5, 2e-5]

visits = [
	locate_x1, locate_x2, locate_x3, locate_x4
]

# bounds on problem
tof_max = 12*365*86400/cparams.tstar
tf = (etf - et0)/cparams.tstar
t0_bnds = [0.5tf, tf]
vinf0_max = 5.0/cparams.vstar

# get objective function
mga1dsm_problem!, lx, ux, ng, lg, ug = Lambert.construct_mga1dsm_problem(
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
visit_nodes, kepler_res, lamb_res, planets = Lambert.view_mga1dsm_problem(
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
	res = Lambert.propagate_arc(lamb)
	plot!(ptraj, res[1,:], res[2,:], c=:orange)
end

# inital and final orbit
cs = [:blue, :crimson, :purple, :brown]
for (i,planet) in enumerate(planets)
	scatter!(ptraj, [planet[1,1]], [planet[2,1]], marker=:circle, c=cs[i])
	plot!(ptraj, planet[1,:], planet[2,:], linestyle=:dash, c=cs[i])
end
# display plot
display(ptraj)

println("Done!")

    