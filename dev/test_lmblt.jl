"""
Test for Lambert-LT problem
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
x02 = [1.5, 0.0, 0.0, 0.0, sqrt(1/1.5), 0.0]
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
n = 3

t0min = 0.0
t0max = psyn
tofmax = period1

lamblt_problem!, lx, ux, ng, lg, ug = construct_lamblt_problem(
	visits, 
	mu,
	n,
	t0min,
	t0max,
	tofmax
)

## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 1500,   # 1500 ~ 2500
    "print_level" => 5,
    "tol" => 1e-4
)

# random initial guess
x0 = random_initial_guess(lx,ux)
println("x0: ", length(x0))

xopt, fopt, info = minimize(lamblt_problem!, x0, ng; 
	lx=lx, 
	ux=ux, 
	lg=lg, 
	ug=ug, 
	solver="ipopt",
	options=ip_options
)

# view results
arcs, dvs, rs = view_lamblt_problem(
	visits, 
	mu,
	n,
	t0min,
	t0max,
	tofmax,
	xopt
)

# propagate orbits for plotting
prop_orb1 = get_plot_orbit(x01, mu)
prop_orb2 = get_plot_orbit(x02, mu)

# key epochs
et0 = xopt[1]
et1 = xopt[1] + xopt[2]

posA_et0, _ = locate_x1(et0)
posA_et1, _ = locate_x1(et1)

posB_et0, _ = locate_x2(et0)
posB_et1, _ = locate_x2(et1)

# prepare plot
ptraj = plot(
	set_aspect=:equal, 
	legend=false, 
	size=(400,400),
	frame_style=:box, 
	ylim=[-3.0,3.0],
	xlim=[-3.0,3.0],
)

cs = palette([:blue, :red], length(rs))

# planets
plot!(ptraj, prop_orb1[1,:], prop_orb1[2,:], c=:black)
plot!(ptraj, prop_orb2[1,:], prop_orb2[2,:], c=:black)

# trajectory
for (i,arc) in enumerate(arcs)
	sol = propagate_arc(arc)
	plot!(ptraj, sol[1,:], sol[2,:], c=:purple)
	# if i > 1
	# 	scatter!(ptraj, [arc.r1[1]], [arc.r1[2]], marker=:circle, c=cs[i])
	# end
end

for (i,r) in enumerate(rs)
	scatter!(ptraj, [r[1]], [r[2]], marker=:circle, c=cs[i])
end
# plot!(ptraj, prop1[1,:], prop1[2,:], c=:blue)
# plot!(ptraj, prop2[1,:], prop2[2,:], c=:red)

display(ptraj)


println("Done!")




