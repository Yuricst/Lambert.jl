"""
Impulsive rendez-vous problems with 2 impulses
"""


"""
	function construct_rdv2imp_problem(
		visits::Vector,
	    mu::Float64=1.0,
	    rev_max::Int=0,
	)

Construct two-impulse rendez-vous problem
"""
function construct_rdv2imp_problem(
	visits::Vector,
    mu::Float64=1.0,
    rev_max::Int=0,
)
	# constraints
	ng = 1   # place-holder constraint
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
	function construct_rdv2imp_problem_tfcon(
		visits::Vector,
		tf_max::Float64,  # final arrival epoch
	    mu::Float64=1.0,
	    rev_max::Int=0,
	)

Construct two-impulse rendez-vous problem with max arrival epoch
"""
function construct_rdv2imp_problem_tfcon(
	visits::Vector,
	tf_max::Float64,  # final arrival epoch
    mu::Float64=1.0,
    rev_max::Int=0,
)
	# constraints
	ng = 1   # FIXME
	lg = [0.0]
	ug = [tf_max]

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
		g[1] = t0+tof
		return dv1+dv2
	end
	return rdv2imp!, ng, lg, ug
end


"""
	function view_rdv2imp_problem(
		x::Vector,
		visits::Vector,
		et0::Float64=0.0,
		tstar::Float64=1.0,
	    vstar::Float64=1.0,
	    use_spice::Bool=true,
	    mu::Float64=1.0,
	    rev_max::Int=0,
	)

View solution of rdv2imp problem
"""
function view_rdv2imp_problem(
	x::Vector,
	visits::Vector,
	et0::Float64=0.0,
	tstar::Float64=1.0,
    vstar::Float64=1.0,
    use_spice::Bool=true,
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

	if use_spice == true
		@printf("Departure epoch:     %s\n", et2utc(et0+t0*tstar,"C",0))
		@printf("Arrival epoch:       %s\n", et2utc(et0+(t0+tof)*tstar,"C",0))
	else
		@printf("Departure epoch:     %1.6e\n",(t0+t0_min)*tstar)
		@printf("Arrival epoch:       %1.6e\n",t0+tof)
	end
	@printf("Time of flight:      %1.6e\n",tof)
	@printf("Initial phase angle: % 3.3f [deg]\n", ϕ0*180/π)
	@printf("Final phase angle:   % 3.3f [deg]\n", ϕf*180/π)
	@printf("Initial ΔV:          %1.6f\n", dv1*vstar)
	@printf("Final ΔV:            %1.6f\n", dv2*vstar)
	@printf("Total ΔV:            %1.6f\n", (dv1+dv2)*vstar)

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


"""
	construct_rdv2imp_spice_problem(
		bodies::Vector{Int}=[3,4],
		et0_str::String="2022-01-01T00:00:00.00", 
		etf_str::String="2024-01-01T00:00:00.00",
		tof_max_day::Float64=365.0,
		canonical_params=nothing,
	)

SPICE-based rendez-vous problem wrapper
"""
function construct_rdv2imp_spice_problem(
	bodies::Vector{Int}=[3,4],
	et0_str::String="2022-01-01T00:00:00.00", 
	etf_str::String="2025-01-01T00:00:00.00",
	tof_max_day::Float64=365.0,
	canonical_params=nothing,
)
	# convert epochs
	et0 = utc2et(et0_str)
	etf = utc2et(etf_str)
	return construct_rdv2imp_spice_problem(
		bodies,
		et0,
		etf,
		tof_max_day,
		canonical_params,
	)
end


"""
	function construct_rdv2imp_spice_problem(
		bodies::Vector{Int},
		et0::Float64, 
		etf::Float64,
		tof_max_day::Float64=365.0,
		canonical_params=nothing,
	)

SPICE-based rendez-vous problem wrapper
"""
function construct_rdv2imp_spice_problem(
	bodies::Vector{Int},
	et0::Float64, 
	etf::Float64,
	tof_max_day::Float64=365.0,
	canonical_params=nothing,
)

	# if canonical parameters is not provided, create Sun-Earth based ones
	if isnothing(canonical_params)
		# dynamics constants
		MU_SUN = 1.3271244004193938e11
		mu =  1.0
		lstar = 1.495978707e8  # 1AU in km
		vstar = sqrt(MU_SUN / lstar)
		tstar = lstar / vstar
		canonical_params = SunEarthCanonical(
			lstar, vstar, tstar, et0
		)
	end

	# scaled time-windows
	t0 = 0.0
	tf = (etf - et0)/canonical_params.tstar

	# construct functions to locate positions
	function locate_x1(epoch::Float64)
		sv = spkssb(bodies[1], et0 + epoch*canonical_params.tstar, "ECLIPJ2000")
		return sv[1:3]/lstar, sv[4:6]/canonical_params.vstar
	end

	function locate_x2(epoch::Float64)
		sv = spkssb(bodies[2], et0 + epoch*canonical_params.tstar, "ECLIPJ2000")
		return sv[1:3]/lstar, sv[4:6]/canonical_params.vstar
	end

	# bounds on variables
	tof_max = tof_max_day*86400/canonical_params.tstar
	lx = [t0; 0.0]
	ux = [tf; tof_max]

	# get objective function
	visits = [locate_x1, locate_x2]
	rendezvous!, ng, lg, ug = construct_rdv2imp_problem_tfcon(
		visits, 
		tf,
		mu,
	)

	return rendezvous!, lx, ux, ng, lg, ug, visits, canonical_params
end




