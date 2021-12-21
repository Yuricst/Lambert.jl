"""
Test for cycler problem
"""

using joptimise
using LinearAlgebra
using Printf
using Plots
#plotly()

push!(LOAD_PATH, "../")
using Lambert

# constants
MU_SUN = 132_712_440_018
MU_EARTH = 398600.44
MU_MARS = 42848.0
AU = 149.6e6
R_EARTH = 6378.0
R_MARS = 3396.0

# helper functions
function get_plot_orbit(state, mu, tf=nothing, steps=500)
	if isnothing(tf) == true
		tf = get_period(state, mu)
	end
	ts = LinRange(0.0, tf, steps)
	prop = zeros(6,steps)
	prop[:,1] = state

	for k = 1:steps-1
		prop[:,k+1] = keplerder_nostm(mu, state, 0.0, ts[k+1], 1.e-14, 20)
	end
	return prop
end


# initial and final condition
x01 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
x02 = [1.5, 0.0, 0.0002, 0.0, sqrt(1/1.5), 0.0]
mu =  1.0
m = 0

period1 = get_period(x01, mu)
period2 = get_period(x02, mu)
psyn = 1/abs(1/period1 - 1/period2)
println("period1: $period1")
println("Synodic period: $psyn")

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

# construct optimization problem
visits = [locate_x1, locate_x2]
cycle_period = psyn * 2
trange = []
body_mus = [MU_EARTH/MU_SUN, MU_MARS/MU_SUN]
min_radii = [1.05*R_EARTH/AU, 1.05*R_MARS/AU]

cycler2_problem!, ng, lg, ug = construct_cycler2_problem(
	visits, 
	cycle_period,
	body_mus,
	min_radii,
	mu,
)

# initial guess
et0min = 0.0
et0max = max(psyn, cycle_period)
x0 = rand(2)
#[2.0764719984822677, 0.4324772298617246]#  # [(et0min+et0max)/2, 0.4]
# bounds on variables
lx = [et0min; 0.0]
ux = [et0max; 1.0]

## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 2500,   # 1500 ~ 2500
    "print_level" => 4,
    "tol" => 1e-6
)

xopt, fopt, info = minimize(cycler2_problem!, x0, ng; 
	lx=lx, 
	ux=ux, 
	lg=lg, 
	ug=ug, 
	solver="ipopt",
	options=ip_options
)

println(xopt)
println(info)

# view results
println("min_radii: $min_radii")
prop1, prop2 = view_cycler2_problem(
	visits, 
	cycle_period,
	body_mus,
	min_radii,
	mu,
	xopt
)

# propagate orbits for plotting
prop_orb1 = get_plot_orbit(x01, mu, cycle_period)
prop_orb2 = get_plot_orbit(x02, mu)

posA1,_ = locate_x1(xopt[1])
posA2,_ = locate_x1(xopt[1] + cycle_period)
posB,_  = locate_x2(xopt[1] + xopt[2]*cycle_period)

# prepare plot
ptraj = plot(
	set_aspect=:equal, 
	legend=false, 
	size=(400,400),
	frame_style=:box, 
	ylim=[-3.0,3.0],
	xlim=[-3.0,3.0],
)
#plot!([-2.0,2.0],[-2.0,2.0],linealpha=0.1)

# planets
plot!(ptraj, prop_orb1[1,:], prop_orb1[2,:], c=:black)
plot!(ptraj, prop_orb2[1,:], prop_orb2[2,:], c=:black)
scatter!(ptraj, [posA1[1]], [posA1[2]], c=:blue, marker=:circle)
scatter!(ptraj, [posA2[1]], [posA2[2]], c=:green, marker=:circle)
scatter!(ptraj, [posB[1]], [posB[2]], c=:red, marker=:circle)

# trajectory
plot!(ptraj, prop1[1,:], prop1[2,:], c=:blue)
plot!(ptraj, prop2[1,:], prop2[2,:], c=:green)

display(ptraj)
println("Done!")




