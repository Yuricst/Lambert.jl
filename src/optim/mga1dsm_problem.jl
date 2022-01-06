"""
MGA problem, no DSM
"""


"""
    unpack_mga1dsm_x(x::Vector, np::Int)

Unpack decision vector for mga1dsm problem
"""
function unpack_mga1dsm_x(x::Vector, np::Int)
    # launch parameters
    t0 = x[1]
    # phase parameters
    p_phases = []
    tofs = []
    for i = 1:np
        if i == 1
            push!(p_phases, x[2+(i-1)*5:1+5i])
            push!(tofs, x[6])
        else
            push!(p_phases, x[7+(i-2)*4:4+6(i-1)])
            push!(tofs, x[6+4(i-1)])
        end
    end
    return t0, p_phases, tofs
end



"""
    function flyby_hyperbolic(
        v0::Vector{Float64}, 
        vp::Vector{Float64},
        rp::Float64,
        β::Float64,
        μ_body::Float64,
        Δv::Float64=0.0,
    )

Patched-conics hyperbolic fly-by
"""
function flyby_hyperbolic(
    v0::Vector{Float64}, 
    vp::Vector{Float64},
    rp::Float64,
    β::Float64,
    μ_body::Float64,
    Δv::Float64=0.0,
)
    vinf0 = v0 - vp
    vinf0_mag = norm(vinf0)
    ecc = 1 + rp/μ_body * vinf0_mag^2
    δ = asin_safe(1/ecc)
    vinf1 = (vinf0_mag+Δv)*[cos(δ), cos(β)*sin(δ), sin(β)*sin(δ)]
    return vinf1 + vp
end


"""
    function construct_mga1dsm_problem(
        visits::Vector, 
        body_mus::Vector,
        body_radii::Vector,
        h_safes::Vector,
        mu::Float64, 
        vinf0_max::Float64=0.0, 
        rev_max::Int=0
    )

Construct MGA-1DSM problem, direct encoding
See for reference: 
- https://esa.github.io/pykep/documentation/trajopt.html#pykep.trajopt.mga_1dsm
- https://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-RPR-MAD-2010-(CambridgePress)GlobalOptimizationAndSpacePruningForSpacecraftTrajectoryDesign.pdf
"""
function construct_mga1dsm_problem(
    visits::Vector, 
    mu::Float64,
    t0_bnds::Vector{Float64},
    tof_max::Float64, 
    vinf0_max::Float64,
    r_bodies::Vector{Float64}=[1.e-6],
    μ_bodies::Vector{Float64}=[1.e-6],
    add_vinf_dep::Bool=false,
    add_vinf_arr::Bool=true,
    rev_max::Int=0,
    η_lb::Float64=0.1,
    η_ub::Float64=0.9,
    tof_min_ratio::Float64=0.005,
)
    # number of phases
    np = length(visits) - 1

    # define minimum tof
    tof_min = tof_min_ratio*tof_max

    # bounds on decision vector
    lx = vcat([t0_bnds[1]], [0.0, 0.0, 0.0, η_lb, tof_min])[:]
    ux = vcat([t0_bnds[2]], [1.0, 1.0, vinf0_max, η_ub, tof_max])[:]
    if np > 1  # append to bounds vector
        for j = 1:np-1
            lx = vcat(lx, [-π, 0.0, η_lb, tof_min])[:]
            ux = vcat(ux, [ π, Inf, η_ub, tof_max])[:]
        end
    end

    # constraints information (on total time of flight)
    ng = 1
    lg = [-Inf]
    ug = [0.0]

    # storage for previous lambert velocity
    lamb_prev_vf = [0.0,0.0,0.0]

    """MGA problem objective function"""
    function mga1dsm_problem!(g, x)
        # unpack decision vector
        t0, p_phases, tofs = unpack_mga1dsm_x(x, np)

        # get nodes to visit
        visit_nodes = []
        for (k, visit) in enumerate(visits)
            if k == 1
                epoch = t0
            else
                epoch = t0 + sum(tofs[1:k-1])
            end
            push!(visit_nodes, visit(epoch))
        end

        # solve consecutive Lambert problems
        if add_vinf_dep == true
            dv_total = 0.0 + p_phases[1][3]  # add v-inf at departure a priori
        else
            dv_total = 0.0
        end
        for k = 1:np
            # unpack phase parameters
            if k == 1
                u, v, vinf0, η, tof = p_phases[k]
            else
                β, rp_rV, η, tof = p_phases[k]
            end

            # get planet positions
            r1vec, v1vec = visit_nodes[k]
            r2vec, v2vec = visit_nodes[k+1]

            # if initial phase, apply v-infinity
            if k == 1
                # convert sphere point picking parameters u and v
                θ = 2π*u
                ϕ = acos_safe(2v - 1) - π/2
                v1vec = v1vec + vinf0*[cos(θ)*cos(ϕ), sin(θ)*cos(ϕ), sin(ϕ)]
            # else, apply fly-by
            else
                v1vec = flyby_hyperbolic(
                    lamb_prev_vf,
                    v1vec, 
                    rp_rV*r_bodies[k],
                    β,
                    μ_bodies[k],
                )
            end

            # propagate forward initial state
            sv_before_dsm = keplerder_nostm(
                mu,
                vcat(r1vec, v1vec)[:], 
                0.0,
                tof*η,
                1e-14,
                20,
            )

            # solve Lambert problem
            res = lambert_fast(
                sv_before_dsm[1:3], 
                r2vec, 
                tof*(1-η), 
                rev_max, 
                mu
            )
            lamb_prev_vf[:] = res.v2

            # append dsm cost
            dsm_cost = norm(res.v1 - sv_before_dsm[4:6])
            dv_total += dsm_cost

            # append final velocity mismatch cost
            if add_vinf_arr == true || k == np
                final_dv = norm(v2vec - res.v2)
                dv_total += final_dv
            end
        end

        # store constraints on total tof
        g[1] = sum(tofs) - tof_max

        # return objective
        return dv_total
    end

    # return generated callable
    return mga1dsm_problem!, lx, ux, ng, lg, ug
end



function view_mga1dsm_problem(
    x::Vector,
    visits::Vector, 
    mu::Float64,
    t0_bnds::Vector{Float64},
    tof_max::Float64, 
    vinf0_max::Float64,
    r_bodies::Vector{Float64}=[1.e-6],
    μ_bodies::Vector{Float64}=[1.e-6],
    add_vinf_dep::Bool=false,
    add_vinf_arr::Bool=true,
    rev_max::Int=0,
    η_lb::Float64=0.1,
    η_ub::Float64=0.9,
)
    # create dummy canonical parameters
    DummyCP = SunEarthCanonical(1.0, 1.0, 1.0, 0.0)
    view_mga1dsm_problem(
        DummyCP,
        x,
        visits,
        mu,
        t0_bnds,
        tof_max, 
        vinf0_max,
        r_bodies,
        μ_bodies,
        add_vinf_dep,
        add_vinf_arr,
        rev_max,
        η_lb,
        η_ub,
    )
end



"""
View result from mga1dsm problem
"""
function view_mga1dsm_problem(
    cparams::CanonicalParameters,
    x::Vector,
    visits::Vector, 
    mu::Float64,
    t0_bnds::Vector{Float64},
    tof_max::Float64, 
    vinf0_max::Float64,
    r_bodies::Vector{Float64}=[1.e-6],
    μ_bodies::Vector{Float64}=[1.e-6],
    add_vinf_dep::Bool=false,
    add_vinf_arr::Bool=true,
    rev_max::Int=0,
    η_lb::Float64=0.1,
    η_ub::Float64=0.9,
)
    # number of phases
    np = length(visits) - 1

    # unpack decision vector
    t0, p_phases, tofs = unpack_mga1dsm_x(x, np)

    # get nodes to visit
    visit_nodes = []
    for (k, visit) in enumerate(visits)
        if k == 1
            epoch = t0
        else
            epoch = t0 + sum(tofs[1:k-1])
        end
        push!(visit_nodes, visit(epoch))
    end

    # storage for analysis
    kepler_res = []
    lamb_res = []
    dsm_vectors = []
    planets = []

    # storage for previous lambert velocity
    lamb_prev_vf = [0.0,0.0,0.0]

    println("***** MGA1DSM transfer *****")
    println("add_vinf_dep: $add_vinf_dep, add_vinf_arr: $add_vinf_arr")
    @printf("Departure epoch:      %s\n", et2utc(cparams.et0+t0*cparams.tstar,"C",0))
    @printf("Total time of flight: %1.6e\n", sum(tofs)*cparams.tstar/86400)
    @printf("Departure ΔV:         %1.6f\n",p_phases[1][3]*cparams.vstar)

    # solve consecutive Lambert problems
    if add_vinf_dep == true
        dv_total = 0.0 + p_phases[1][3]  # add v-inf at departure a priori
    else
        dv_total = 0.0
    end
    for k = 1:np
        # unpack phase parameters
        if k == 1
            u, v, vinf0, η, tof = p_phases[k]
        else
            β, rp_rV, η, tof = p_phases[k]
        end

        # get planet positions
        r1vec, v1vec = visit_nodes[k]
        r2vec, v2vec = visit_nodes[k+1]

        # if initial phase, apply v-infinity
        if k == 1
            # convert sphere point picking parameters u and v
            θ = 2π*u
            ϕ = acos_safe(2v - 1) - π/2
            v1vec = v1vec + vinf0*[cos(θ)*cos(ϕ), sin(θ)*cos(ϕ), sin(ϕ)]
        # else, apply fly-by
        else
            v1vec = flyby_hyperbolic(
                lamb_prev_vf,
                v1vec, 
                rp_rV*r_bodies[k],
                β,
                μ_bodies[k],
            )
        end

        # propagate forward initial state
        sv_before_dsm = keplerder_nostm(
            mu,
            vcat(r1vec, v1vec)[:], 
            0.0,
            tof*η,
            1e-14,
            20,
        )
        # propagate over interval for storing and plotting
        kepprop = keplerder_nostm(
            mu, 
            vcat(r1vec, v1vec)[:], 
            0.0, 
            LinRange(0.0, tof*η, 500), 
            1.e-14, 
            20
        )

        # solve Lambert problem
        res = lambert_fast(
            sv_before_dsm[1:3], 
            r2vec, 
            tof*(1-η), 
            rev_max, 
            mu
        )
        lamb_prev_vf[:] = res.v2

        # append cost
        dsm_cost = norm(res.v1 - sv_before_dsm[4:6])
        dv_total += dsm_cost

        # append final velocity mismatch cost
        if add_vinf_arr == true && k == np
            final_dv = norm(v2vec - res.v2)
            dv_total += final_dv
        end

        # append to lists
        push!(kepler_res, kepprop)
        push!(lamb_res, res)
        push!(dsm_vectors, res.v1 - sv_before_dsm[4:6])

        # propagate initial and final orbit of phase
        r1vec, v1vec = visit_nodes[k]  # re-write for planet plotting
        prop_r1 = keplerder_nostm(
            mu,
            vcat(r1vec, v1vec)[:],
            0.0,
            LinRange(0.0, get_period(vcat(r1vec, v1vec)[:], mu), 500),
        )
        push!(planets, prop_r1)
        if k == np
            prop_r2 = keplerder_nostm(
                mu,
                vcat(r2vec, v2vec)[:],
                0.0,
                LinRange(0.0, get_period(vcat(r2vec, v2vec)[:], mu), 500),
            )
            push!(planets, prop_r2)
        end

        @printf("Leg %1.0f: \n",k)
        @printf("   DSM cost:       %1.6e\n",dsm_cost*cparams.vstar)
        @printf("   Time of flight: %1.6e\n",tof*cparams.tstar/86400)
        @printf("   Time until DSM: %1.6e\n",tof*η*cparams.tstar/86400)
        @printf("   η:              %1.6e\n",η)
        if add_vinf_arr == true && k == np
            @printf("Arrival ΔV:           %1.6f\n",final_dv*cparams.vstar)
        end
    end
    @printf("Total ΔV:             %1.6f\n",dv_total*cparams.vstar)
    # return results for plotting etc
    return visit_nodes, kepler_res, lamb_res, planets
end

