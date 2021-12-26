"""Functions for lambert-arc based low-thrust approx."""


"""

Construct optimization problem for lambert-arc based low-thrust transfers
Encoding of decision vector: 
	`x = [t0, tof] + np*n*[k,ν,τ]`
"""
function construct_lamblt_problem(
	visits::Vector, 
	mu::Float64,
	n::Int,
	t0min,
	t0max,
	tofmax,
)
	# number of phases
	np = length(visits)-1

	# get periods of each body
	periods = []
	smas = []
	for visit in visits
		rtmp, vtmp = visit(t0min)
		push!(periods, get_period(vcat(rtmp,vtmp)[:], mu))
		push!(smas, get_semiMajorAxis(vcat(rtmp,vtmp)[:], mu))
	end
	
	"""Lambert-arc low-thrust problem objective function"""
	function lamblt_problem!(g, x)
		# unpack decision vector
		t0, tof = x[1:2]

		# list of delta-Vs
		dvs = []

		# iterate for each phase
		for ip = 1:np
			# epochs of nodes
			t1 = t0
			t2 = t0 + tof

			r1, v1 = visits[ip](t1)
			r2, v2 = visits[ip+1](t2)
			# convert to spherical coordinates
			r1_sph = pos_cartesian2spherical(r1)
			r2_sph = pos_cartesian2spherical(r2)
			# discretize space by n
			θ1, θ2 = r1_sph[2], r2_sph[2]
			if θ1 > θ2
				θ2 = θ2 + 2π
			end
			θs = LinRange(θ1, θ2, n+2)
			ϕs = LinRange(r1_sph[3], r2_sph[3], n+2)

			# re-construct intermediate position vectors
			rs = [r1]
			for ii = 1:n
				# extract info from decision vector
				k = x[2 + 3*ip*ii-2]
				ν = x[2 + 3*ip*ii-1]
				# re-construct intermediate position vectors
				r_inter_sphe = [smas[ip]*k, θs[ii+1], ϕs[ii+1]]
				r_inter_sphe_mod = rotmat2(ν) * r_inter_sphe
				r_inter_cart = pos_spherical2cartesian(r_inter_sphe_mod)
				push!(rs, r_inter_cart)
			end
			push!(rs, r2)

			# solve Lambert problem consecutively
			v_prev = v1
			tof_ramining = tof
			for ii = 1:n+1
				# extract info from decision vector
				if ii < n+1
					τ = x[2 + 3*ip*ii]
					dt = tof_ramining*τ
					tof_ramining = tof_ramining - dt
				else
					dt = tof_ramining
				end

				if dt < periods[ip]
					m = 0
				else
					m = 1
				end

				arc = lambert_fast(rs[ii], rs[ii+1], dt, m, mu)
				push!(dvs, norm(arc.v1 - v_prev))
				# update previous delta-V
				v_prev = arc.v2
			end
		end
		# place holder on constraint
		g[1] = 0.0
		# return sum of delta-V as minimization objective
		return sum(dvs)
	end
	# constraints bounds
	lg = [0.0]
	ug = [0.0]
	ng = 1
	# decision variables bounds
	lx = [t0min, 0.0]
	ux = [t0max, tofmax] 
	for idx = 1:np*n
		push!(lx, 0.5)  # lower-bound on k
		push!(ux, 2.0)   # upper-bound on k
		push!(lx, -π/2)  # lower-bound on ν
		push!(ux,  π/2)  # upper-bound on ν
		push!(lx, 0.0)   # lower-bound on τ
		push!(ux, 1.0)   # upper-bound on τ
	end
	return lamblt_problem!, lx, ux, ng, lg, ug
end


function view_lamblt_problem(
	visits::Vector, 
	mu::Float64,
	n::Int,
	t0min,
	t0max,
	tofmax,
	x::Vector,
)	
	# number of phases
	np = length(visits)-1

	# get periods of each body
	periods = []
	smas = []
	for visit in visits
		rtmp, vtmp = visit(t0min)
		push!(periods, get_period(vcat(rtmp,vtmp)[:], mu))
		push!(smas, get_semiMajorAxis(vcat(rtmp,vtmp)[:], mu))
	end
	
	# unpack decision vector
	t0, tof = x[1:2]

	# print information
	@printf("Departure:      %1.4e\n", t0)
	@printf("Time of flight: %1.4e\n", tof)
	tof_ramining = tof
	for ip = 1:np
		for ii = 1:n+1
			if ii < n+1
				k = x[2 + 3*ip*ii-2]
				ν = x[2 + 3*ip*ii-1]
				τ = x[2 + 3*ip*ii]
				dt = tof_ramining*τ
				tof_ramining = tof_ramining - dt
			else
				dt = tof_ramining
			end
			if dt < periods[ip]
				m = 0
			else
				m = 1
			end

			@printf("Segment %2.0f:\n",ii)
			if ii < n+1
				@printf("   k:              %1.4e\n",k)
				@printf("   ν [deg]:        %3.4f\n",ν*180/π)
			end
			@printf("   Time of flight: %1.4e\n", dt)
			@printf("   Lambert rev:    %1.0f\n", m)
		end
	end

	# initialize list of delta-Vs
	arcs = []
	dvs = []
	rs = []

	# iterate for each phase
	for ip = 1:np
		# epochs of nodes
		t1 = t0
		t2 = t0 + tof

		r1, v1 = visits[ip](t1)
		r2, v2 = visits[ip+1](t2)
		# convert to spherical coordinates
		r1_sph = pos_cartesian2spherical(r1)
		r2_sph = pos_cartesian2spherical(r2)
		# discretize space by n
		θ1, θ2 = r1_sph[2], r2_sph[2]
		if θ1 > θ2
			θ2 = θ2 + 2π
		end
		θs = LinRange(θ1, θ2, n+2)
		ϕs = LinRange(r1_sph[3], r2_sph[3], n+2)

		# re-construct intermediate position vectors
		push!(rs, r1)
		for ii = 1:n
			# extract info from decision vector
			k = x[2 + 3*ip*ii-2]
			ν = x[2 + 3*ip*ii-1]
			# re-construct intermediate position vectors
			r_inter_sphe = [smas[ip]*k, θs[ii+1], ϕs[ii+1]]
			r_inter_sphe_mod = rotmat2(ν) * r_inter_sphe
			r_inter_cart = pos_spherical2cartesian(r_inter_sphe_mod)
			push!(rs, r_inter_cart)
		end
		push!(rs, r2)

		# solve Lambert problem consecutively
		v_prev = v1
		tof_ramining = tof
		for ii = 1:n+1
			# extract info from decision vector
			if ii < n+1
				τ = x[2 + 3*ip*ii]
				dt = tof_ramining*τ
				tof_ramining = tof_ramining - dt
			else
				dt = tof_ramining
			end

			if dt < periods[ip]
				m = 0
			else
				m = 1
			end

			arc = lambert_fast(rs[ii], rs[ii+1], dt, m, mu)
			push!(dvs, norm(arc.v1 - v_prev))
			push!(arcs, arc)
			# update previous delta-V
			v_prev = arc.v2
		end
	end
	# return arcs and dvs
	return arcs, dvs, rs
end


"""
Convert cartesian position vector to spherical
"""
function pos_cartesian2spherical(r_cartesian, to2Pi::Bool=true)
	# unpack position-vector
    x,y,z = r_cartesian
    r = sqrt(x^2 + y^2 + z^2)
    θ = atan(y,x)
    if to2Pi == true && θ < 0.0
    	θ = 2π - θ
    end
    r_spherical = [
        r,
        atan(y,x),
        asin(z/r),
    ]
    return r_spherical
end


"""
Convert spherical position vector to cartesian
"""
function pos_spherical2cartesian(r_spherical)
	# unpack state-vector
    ρ, θ, ϕ = r_spherical
    r_cartesian = [
        ρ*cos(ϕ)*cos(θ),
        ρ*cos(ϕ)*sin(θ),
        ρ*sin(ϕ)
    ]
    return r_cartesian
end


"""
Rotation matrix about second axis
"""
function rotmat2(ϕ)
	c = cos(ϕ)
	s = sin(ϕ)
	return [c 0 s; 0 1 0; -s 0 c]
end


