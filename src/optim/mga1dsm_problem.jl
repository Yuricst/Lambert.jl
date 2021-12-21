"""
MGA problem, no DSM
"""


function unpack_mga1dsm_x(x::Vector, np::Int)
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
Construct MGA-1DSM problem
"""
function construct_mga1dsm_problem(visits::Vector, mu, vinf0_max, rev_max = 1)

    # number of phases
    np = length(visits) - 1

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
        dv_total = 0.0
        for k = 1:np
            # get position vectors
            r1vec = visit_nodes[k][1]
            r2vec = visit_nodes[k][1]

            # solve Lambert problem
            res_list = []
            for klmb = 1:rev_max+1
                push!(res_list, lambert_fast(r1vec, r2vec, tofs[k], klmb - 1, mu))
            end

            # compute cost
            dv_scen = []
            for res in res_list
                # for launch node, simply compute launch v-infinity
                if k == 1
                    v1vec = visit_nodes[k][2]
                    push!(dv_scen, res.v1 - v1vec)

                    # for fly-by node, use fly-by
                else
                end
            end

            # keep minimum
            dv_total += minimum(dv_scen)
        end

        # store constraints

        # return objective
        return
    end

    # return generated callable
    return mga_problem!
end
