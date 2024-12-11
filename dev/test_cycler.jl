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

# bounds on variables
et0min = 0.0
et0max = max(psyn, cycle_period)
# bounds on variables
lx = [et0min; 0.0]
ux = [et0max; 1.0]

## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 500,   # 1500 ~ 2500
    "print_level" => 0,
    "tol" => 1e-6
)

# number of retry
retry_num = 10
# initialize storage
sols = []
for i = 1:retry_num
	# random initial guess
	x0 = rand(2)
	#[2.0764719984822677, 0.4324772298617246]#  # [(et0min+et0max)/2, 0.4]

	xopt, fopt, info = minimize(cycler2_problem!, x0, ng; 
		lx=lx, 
		ux=ux, 
		lg=lg, 
		ug=ug, 
		solver="ipopt",
		options=ip_options
	)

	if info == :Solve_Succeeded
		println("Success at try $i")
		push!(sols, [xopt, fopt, info])
		break
	end
end

if length(sols)>0
	# if successful, extract result
	xopt, fopt, info = sols[1]

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

	# key epochs
	et0 = xopt[1]
	et1 = xopt[1] + cycle_period*xopt[2]
	et2 = xopt[1] + cycle_period

	posA_et0, _ = locate_x1(et0)
	posA_et1, _ = locate_x1(et1)
	posA_et2, _ = locate_x1(et2)

	posB_et0, _ = locate_x2(et0)
	posB_et1, _ = locate_x2(et1)
	posB_et2, _ = locate_x2(et2)

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

	cs = palette([:purple, :green], 3)

	# planets
	plot!(ptraj, prop_orb1[1,:], prop_orb1[2,:], c=:black)
	plot!(ptraj, prop_orb2[1,:], prop_orb2[2,:], c=:black)

	# planets
	scatter!(ptraj, [posA_et0[1]], [posA_et0[2]], c=cs[1], marker=:circle)
	scatter!(ptraj, [posA_et1[1]], [posA_et1[2]], c=cs[2], marker=:circle)
	scatter!(ptraj, [posA_et2[1]], [posA_et2[2]], c=cs[3], marker=:circle)

	scatter!(ptraj, [posB_et0[1]], [posB_et0[2]], c=cs[1], marker=:circle)
	scatter!(ptraj, [posB_et1[1]], [posB_et1[2]], c=cs[2], marker=:circle)
	scatter!(ptraj, [posB_et2[1]], [posB_et2[2]], c=cs[3], marker=:circle)

	# initial and final phase
	plot!(ptraj, [0, posA_et0[1]], [0, posA_et0[2]], c=cs[1], linewidht=0.5)
	plot!(ptraj, [0, posB_et0[1]], [0, posB_et0[2]], c=cs[1], linewidht=0.5)

	plot!(ptraj, [0, posA_et2[1]], [0, posA_et2[2]], c=cs[3], linewidht=0.5)
	plot!(ptraj, [0, posB_et2[1]], [0, posB_et2[2]], c=cs[3], linewidht=0.5)

	# trajectory
	plot!(ptraj, prop1[1,:], prop1[2,:], c=:blue)
	plot!(ptraj, prop2[1,:], prop2[2,:], c=:red)

	display(ptraj)
end

println("Done!")




