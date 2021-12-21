"""
MGA problem, no DSM
"""


function unpack_mga_x(x::Vector, np::Int)
    # launch parameters
    t0 = x[1]
    # phase parameters
    tofs = []
    pumps = []
    cranks = []
    for i = 1:np
        push!(tofs, x[1+i+3*(i-1)])
        if i < np   # if not final node, also store pump and crank angles
            push!(pumps, x[1+i+1+3*(i-1)])
            push!(cranks, x[1+i+2+3*(i-1)])
        end
    end
    return t0, tofs, pumps, cranks
end


"""
Construct MGA problem
"""
function construct_mga_problem(
    visits::Vector,
    body_mus::Vector,
    body_radii::Vector,
    h_safes::Vector,
    mu::Float64,
    rev_max = 1,
)

    # number of phases
    np = length(visits) - 1
    # number of constraints
    ng = max(0, np - 1)

    """MGA problem objective function"""
    function mga_problem!(g, x)
        # unpack decision vector
        t0, tofs, pumps, cranks = unpack_mga_x(x, np)

        # get states to visit
        visit_nodes = []
        for (k, visit) in enumerate(visits)
            push!([visit(t0 + sum(tofs[1:k]))])
        end

        # solve consecutive Lambert problems
        dv_total = 0.0  # initialize storage
        vinf_prev = [0.0, 0.0, 0.0]  # initialize storage
        for k = 1:np
            # get position vectors
            r1vec = visit_nodes[k][1]
            r2vec = visit_nodes[k][1]

            # solve Lambert problem
            res = lambert_fast(r1vec, r2vec, tofs[k], 0, mu)

            if k == 1
                dv_total = res.v1 - v1vec
            end

            # res_list = []
            # for klmb = 1:rev_max+1
            # 	push!(res_list, lambert_fast(r1vec, r2vec, tofs[k], klmb-1, mu))
            # end

            # # compute cost from launch node
            # dv_scen = []
            # if k == 1
            # 	for res in res_list
            # 		v1vec = visit_nodes[k][2]
            # 		push!(dv_scen, res.v1 - v1vec)
            # 	end
            # end
            # # keep minimum
            # dv_total += minimum(dv_scen)

            # compute fly-by constraint
            if k > 1
                # out-going vinf from r1
                vinf_out = res.v1 - v1vec
                # turn-angle
                δ = acos_safe(dot(vinf_in, vinf_out) / (norm(vinf_in) * norm(vinf_out)))
                # fly-by constraint
                g[k-1] =
                    (body_radii[k+1] + h_safes[k+1]) -
                    body_mus[k+1] / (norm(vinf_in))^2 * (1 / sin(δ / 2) - 1)
            end

            # store incoming vinf for next iteration
            vinf_in = res.v2 - v2vec
        end
        # return objective
        return dv_total
    end

    # return generated callable, number of constraints
    return mga_problem!, ng
end
