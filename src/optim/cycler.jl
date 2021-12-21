"""Functions for generating cyclers"""


function unpack_cycler2_x(x::Vector)
	# launch epoch
	et0 = x[1]
	# time of flight fraction from body A -> B
	η = x[2]
end


"""
Construct optimization problem to find 2-body cycler
Encoding of decision vector: 
	`x = [et0, η]`

Constraint:
	- 
	
Objective: 
	`(v∞)^2 at launch`
"""
function construct_cycler2_problem(
	visits::Vector, 
	cycle_period::Float64,
	body_mus::Vector,
	min_radii::Vector,
	mu::Float64,
	m1::Int=0,
	m2::Int=1,
)

	# get minimum fly-by radii
	r1min, r2min = min_radii[1], min_radii[2]

	# number and bounds on constraints
	ng = 4
	lg = [0.0, 0.0, -Inf, -Inf]
	ug = [0.0, 0.0, 0.0, 0.0]

	"""Cycler-2 problem objective function (sequence A-B-A)"""
	function cycler2_problem!(g, x)
		# unpack decision vector
		et0, η = x[1], x[2]
		tof_AB = cycle_period*η
		tof_BA = cycle_period*(1.0-η)

		# solve Lambert problem A -> B
		rA1, vA1 = visits[1](et0)
		rB1, vB1 = visits[2](et0 + tof_AB)
		res1 = lambert_fast(rA1, rB1, tof_AB, m1, mu)

		# solve Lambert problem B -> A
		_, vB2   = visits[2](et0 + tof_AB)
		rA2, vA2 = visits[1](et0 + cycle_period)
		res2 = lambert_fast(rB1, rA2, tof_BA, m2, mu)

		# v-infinity of first leg
		vinf1_A = res1.v1 - vA1
		vinf1_B = res1.v2 - vB1

		# v-infinity of second leg
		vinf2_B = res2.v1 - vB2
		vinf2_A = res2.v2 - vA2

		# compute fly-by angles
		δA = acos(dot(vinf1_A, vinf2_A) / (norm(vinf1_A) * norm(vinf2_A)))
		δB = acos(dot(vinf1_B, vinf2_B) / (norm(vinf1_B) * norm(vinf2_B)))

		# equality constraint on velicty magnitudes
		g[1] = norm(res1.v1) - norm(res2.v2)   # match at A
		g[2] = norm(res1.v2) - norm(res2.v1)   # match at B

		# inequality constraints on fly-by angles
		g[3] = min_radii[1] - body_mus[1]/(norm(vinf1_A))^2*(1/sin(δA/2) - 1)
		g[4] = min_radii[2] - body_mus[2]/(norm(vinf1_B))^2*(1/sin(δB/2) - 1)

        # minimize energy of launch orbit
        return norm(vinf1_A)^2
	end
	return cycler2_problem!, ng, lg, ug
end



function view_cycler2_problem(
	visits::Vector,
	cycle_period::Float64, 
	body_mus::Vector,
	min_radii::Vector,
	mu::Float64,
	x::Vector,
	m1::Int=0,
	m2::Int=1,
	steps::Int=500
	)
	# unpack decision vector
	et0, η = x[1], x[2]
	tof_AB = cycle_period*η
	tof_BA = cycle_period*(1.0-η)

	# solve Lambert problem A -> B
	rA1, vA1 = visits[1](et0)
	rB1, vB1 = visits[2](et0 + tof_AB)
	res1 = lambert_fast(rA1, rB1, tof_AB, m1, mu)

	# solve Lambert problem B -> A
	_, vB2   = visits[2](et0 + tof_AB)
	rA2, vA2 = visits[1](et0 + cycle_period)
	res2 = lambert_fast(rB1, rA2, tof_BA, m2, mu)

	# v-infinity of first leg
	vinf1_A = res1.v1 - vA1
	vinf1_B = res1.v2 - vB1

	# v-infinity of second leg
	vinf2_B = res2.v1 - vB2
	vinf2_A = res2.v2 - vA2

	# compute fly-by angles
	δA = acos(dot(vinf1_A, vinf2_A) / (norm(vinf1_A) * norm(vinf2_A)))
	δB = acos(dot(vinf1_B, vinf2_B) / (norm(vinf1_B) * norm(vinf2_B)))

	# inequality constraints on fly-by angles
	rfbA = body_mus[1]/(norm(vinf1_A))^2*(1/sin(δA/2) - 1)
	rfbB = body_mus[2]/(norm(vinf1_B))^2*(1/sin(δB/2) - 1)

	@printf("δ at A: %3.4f [deg]\n", δA*180/π)
	@printf("δ at B: %3.4f [deg]\n", δB*180/π)
	@printf("Fly-by radius at A: %1.4e\n", rfbA)
	@printf("Fly-by radius at B: %1.4e\n", rfbB)

	# propagate trajectory
	xA1 = vcat(rA1, res1.v1)[:]
	xB2 = vcat(rB1, res2.v1)[:]
	t1s = LinRange(0.0, tof_AB, steps)
	t2s = LinRange(0.0, tof_BA, steps)

	# propagate array
	prop1 = zeros(6,steps)
	prop2 = zeros(6,steps)
	prop1[:,1] = xA1
	prop2[:,1] = xB2

	for k = 1:steps-1
		prop1[:,k+1] = keplerder_nostm(mu, xA1, 0.0, t1s[k+1], 1.e-14, 20)
		prop2[:,k+1] = keplerder_nostm(mu, xB2, 0.0, t2s[k+1], 1.e-14, 20)
	end

	return prop1, prop2
end


