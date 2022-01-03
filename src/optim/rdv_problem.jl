"""
Construct two-impulse rendez-vous problem
"""
function construct_rdv2imp_problem(
	visits::Vector,
    mu::Float64=1.0,
    rev_max::Int=0,
)

	# constraints
	ng = 1   # FIXME
	lg = [0.0]
	ug = [0.0]

	"""Two-impulse rdv problem"""
	function rdv2imp!(g,x)
		# unpack initial vector
		t0, tof = x
		# get initial and final positions
		r1vec, v1vec = visits[1](t0)
		r2vec, v2vec = visits[2](t0+tof)
		# solve Lambert
		res = lambert_fast(r1vec, r2vec, tof, rev_max, mu)
		# compute Delta-V costs
		dv1 = norm(res.v1 - v1vec)
		dv2 = norm(v2vec - res.v2)
		# store constraints (if any)
		g[1] = 0.0
		return dv1+dv2
	end

	return rdv2imp!, ng, lg, ug
end


"""
View solution of rdv2imp problem
"""
function view_rdv2imp_problem(
	x::Vector,
	visits::Vector,
    mu::Float64=1.0,
    rev_max::Int=0,
)
	# unpack initial vector
	t0, tof = x
	# get initial and final positions
	r1vec, v1vec = visits[1](t0)
	r2vec, v2vec = visits[2](t0+tof)
	# solve Lambert
	res = lambert_fast(r1vec, r2vec, tof, rev_max, mu)
	# compute Delta-V costs
	dv1 = norm(res.v1 - v1vec)
	dv2 = norm(v2vec - res.v2)

	# compute phase angles between the planets
	r2_t0, _ = visits[2](t0)
	r1_tf, _ = visits[1](t0+tof)
	ϕ0 = acos_safe(dot(r1vec,r2_t0)/(norm(r1vec)*norm(r2_t0)))
	ϕf = acos_safe(dot(r1_tf,r2vec)/(norm(r1_tf)*norm(r2vec)))

	@printf("Departure epoch:     %1.6e\n",t0)
	@printf("Arrival epoch:       %1.6e\n",t0+tof)
	@printf("Time of flight:      %1.6e\n",tof)
	@printf("Initial phase angle: % 3.3f [deg]\n", ϕ0*180/π)
	@printf("Final phase angle:   % 3.3f [deg]\n", ϕf*180/π)

	# propagate transfer arc
	prop_traj = keplerder_nostm(
		mu,
		vcat(r1vec, res.v1)[:], 
		0.0,
		LinRange(0.0, tof, 500),
	)

	# propagate initial and final orbit
	prop_r1 = keplerder_nostm(
		mu,
		vcat(r1vec, v1vec)[:],
		0.0,
		LinRange(0.0, get_period(vcat(r1vec, v1vec)[:], mu), 500),
	)

	prop_r2 = keplerder_nostm(
		mu,
		vcat(r2vec, v2vec)[:],
		0.0,
		LinRange(0.0, get_period(vcat(r2vec, v2vec)[:], mu), 500),
	)
	return prop_traj, prop_r1, prop_r2
end

