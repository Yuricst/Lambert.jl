"""
Demonstration of simple two-impulse optimization with Lambert solver
"""

using SPICE
using joptimise
using LinearAlgebra
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

furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))

# dynamics constants
MU_SUN = 1.3271244004193938e11
mu =  1.0
lstar = 1.495978707e8  # 1AU in km
vstar = sqrt(MU_SUN / lstar)
tstar = lstar / vstar

# initial and final epoch
et0_str = "2032-01-01T00:00:00.00"
etf_str = "2034-01-01T00:00:00.00"
et0 = utc2et(et0_str)
etf = utc2et(etf_str)

# scaled time-windows
t0 = 0.0
tf = (etf - et0)/tstar

function locate_x1(epoch::Float64)
	# scale epoch
	et = et0 + epoch*tstar
	# get position and velocity vector at epoch
	sv = spkssb(3, et, "ECLIPJ2000")
	return sv[1:3]/lstar, sv[4:6]/vstar
end

function locate_x2(epoch::Float64)
	# scale epoch
	et = et0 + epoch*tstar
	# get position and velocity vector at epoch
	sv = spkssb(4, et, "ECLIPJ2000")
	return sv[1:3]/lstar, sv[4:6]/vstar
end

# initial guess
x0 = rand(2)
# bounds on variables
tof_max = 365*86400/tstar
lx = [t0; 0.0]
ux = [tf; tof_max]

# get objective function
visits = [locate_x1, locate_x2]
rendezvous!, ng, lg, ug = construct_rdv2imp_problem(visits, mu)

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

println(info)

# view result
prop_traj, prop_r1, prop_r2 = view_rdv2imp_problem(
	xopt, visits, et0, tstar, vstar)

# prepare plot
ptraj = plot(
	size=(500,500),
	set_aspect=:equal, 
	legend=false, 
	frame_style=:box, 
)

# intial and final position
r1vec = prop_r1[:,1]
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



