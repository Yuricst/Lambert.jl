"""
Function for propagate using Kepler's equation with Der's formulation
from G. Der (1996).
Yuri Shimane, 2021.03.23
"""


"""Hyperbolic sine function"""
function hypertrig_s(z::Float64)
    if z > 0.0
        s = (sqrt(z) - sin(sqrt(z))) / (sqrt(z))^3
    elseif z < 0.0
        s = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z))^3
    else
        s = 1.0 / 6.0
    end
    return s
end


"""Hyperbolic cosine function"""
function hypertrig_c(z::Float64)
    if z > 0.0
        c = (1.0 - cos(sqrt(z))) / z
    elseif z < 0.0
        c = (cosh(-z) - 1) / (-z)
    else
        c = 0.5
    end
    return c
end


"""Lagrange parameter functions"""
function universal_functions(x::Float64, alpha::Float64)
    """Function computes U0 ~ U3 from G. Der 1996"""
    # evaluate hypertrig function
    S = hypertrig_s(alpha * x^2)
    C = hypertrig_c(alpha * x^2)
    # parameters
    u0 = 1 - alpha * x^2 * C
    u1 = x * (1 - alpha * x^2 * S)
    u2 = x^2 * C
    u3 = x^3 * S
    return u0, u1, u2, u3
end


"""Function computes Kepler's time equation and its derivatives in G. Der 1996 form"""
function kepler_der(
    x::Float64,
    alpha::Float64,
    t::Float64,
    t0::Float64,
    mu::Float64,
    r0::Float64,
    sigma0::Float64,
)
    u0, u1, u2, u3 = universal_functions(x, alpha)
    Fun = r0 * u1 + sigma0 * u2 + u3 - sqrt(mu) * (t - t0)
    dF = r0 * u0 + sigma0 * u1 + u2
    d2F = sigma0 * u0 + (1 - alpha * r0) * u1
    return Fun, dF, d2F
end


"""Function computes Lagrange coefficients (as functionined in G. Der 1996)"""
function lagrange_coefficients_der(
    mu::Float64,
    alpha::Float64,
    r0::Float64,
    v0::Float64,
    sigma0::Float64,
    u0::Float64,
    u1::Float64,
    u2::Float64,
    u3::Float64,
    r::Float64,
    sigma::Float64,
)
    # scalar function
    f = 1.0 - u2 / r0
    g = r0 * u1 / sqrt(mu) + sigma0 * u2 / sqrt(mu)
    fdot = -sqrt(mu) / (r * r0) * u1
    gdot = 1.0 - u2 / r
    return f, g, fdot, gdot
end


"""Scalar coefficients of STM dyads (eqn.18)"""
function der_stm_coefs(
    mu::Float64,
    alpha::Float64,
    r0::Float64,
    v0::Float64,
    r::Float64,
    u0::Float64,
    u1::Float64,
    u2::Float64,
    u3::Float64,
    sigma0::Float64,
    sigma::Float64,
)
    # first row coefficients
    c11 =
        1 / (alpha * r * r0^2) * (3 * u1 * u3 + (alpha * r0 - 2) * u2^2) +
        u1^2 / r +
        u2 / r0
    c12 = v0 * u1 * u2 / (r * sqrt(mu))
    c13 =
        v0 / (alpha * r * r0^2 * sqrt(mu)) * (
            r0 * u1 * u2 + 2 * sigma0 * u2^2 + 3 * u2 * u3 - 3 * r * u3 +
            alpha * r0^2 * u1 * u2
        )
    c14 = v0^2 * u2^2 / (r * mu)
    # second row
    c21 = r0 * u1 * u2 / (r * sqrt(mu))
    c22 = v0 / (alpha * r * mu) * (3 * u1 * u3 + (alpha * r0 - 2) * u2^2)
    c23 = r0 * v0 * u2^2 / (r * mu)
    c24 =
        v0^2 / (alpha * r * mu * sqrt(mu)) *
        (r0 * u1 * u2 + 2 * sigma0 * u2^2 + 3 * u2 * u3 - 3 * r * u3)
    # third row
    c31 =
        sqrt(mu) / (alpha * r^3 * r0^2) * (
            r * (3 * u0 * u3 - u1 * u2 + 2 * alpha * r0 * u1 * u2) -
            sigma * (3 * u1 * u3 - 2 * u2^2 + alpha * r0 * u2^2)
        ) + sqrt(mu) / (r^3 * r0) * (2 * r * r0 * u0 * u1 + r^2 * u1 - sigma * r0 * u1^2)
    c32 = v0 / r^3 * (r * (u0 * u2 + u1^2) - sigma * u1 * u2)
    c33 =
        -3 * v0 * u2 / (alpha * r * r0^2) +
        v0 / (alpha * r^2 * r0^2) *
        (4 * sigma0 * u1 * u2 + r0 * (u0 * u2 + u1^2) - 3 * (u1 * u3 + u2^2)) -
        sigma * v0 * u2 / (alpha * r^3 * r0^2) * (r0 * u1 + 2 * sigma0 * u2 + 3 * u3) +
        v0 / r^3 * (r * (u0 * u2 + u1^2) - sigma * u1 * u2)
    c34 = v0^2 / (r^3 * sqrt(mu)) * (2 * r * u1 * u2 - sigma * u2^2)
    # fourth row
    c41 = r0 / r^3 * (r * (u0 * u2 + u1^2) - sigma * u1 * u2)
    c42 =
        v0 / (alpha * r^3 * sqrt(mu)) * (
            r * (3 * u0 * u3 - u1 * u2 + 2 * alpha * r0 * u1 * u2) -
            sigma * (3 * u1 * u3 - 2 * u2^2 + alpha * r0 * u2^2)
        )
    c43 = r0 * v0 / (r^3 * sqrt(mu)) * (2 * r * u1 * u2 - sigma * u2^2)
    c44 =
        -3 * v0^2 * u2 / (alpha * r * mu) +
        v0^2 / (alpha * r^2 * mu) *
        (4 * sigma0 * u1 * u2 + r0 * (u0 * u2 + u1^2) + 3 * (u1 * u3 + u2^2)) -
        sigma * v0^2 * u2 / (alpha * r^3 * mu) * (r0 * u1 + 2 * sigma0 * u2 + 3 * u3)
    return c11, c12, c13, c14, c21, c22, c23, c24, c31, c32, c33, c34, c41, c42, c43, c44
end


"""
    keplerder(mu::Float64, state0::Vector{Float64}, t0::Float64, t::Float64, tol::Float64=1e-12, maxiter::Int=10)

Function computes position at some future time t, via Keplerian time-propagation with G. Der formulation

# Arguments
    - mu::Float64: gravitational parameter
    - state0 (array): initial state [x,y,z,vx,vy,vz]
    - t0::Float64: initial time at state0
    - t::Float64: final time
    - tol::Float64: tolerance on Laguerre-correction
    - maxiter::Int: max allowed iteration allowed for Laguerre-correction

# Returns
    - (tuple): final state, final STM
"""
function keplerder(
    mu::Float64,
    state0::Vector{Float64},
    t0::Float64,
    t::Float64,
    tol::Float64 = 1e-12,
    maxiter::Int = 10,
)
    # ------------------------------------------ #
    # SET-UP PROBLEM
    r0, v0 = state0[1:3], state0[4:6]
    sma = get_semiMajorAxis(state0, mu)
    alpha = 1.0 / sma
    sigma0 = dot(r0, v0) / sqrt(mu)

    # ------------------------------------------ #
    # ITERATION WITH LAGUERRE-CONWAY METHOD
    # initial guess based on eccentricity
    ecc = norm(GALT.get_eccentricity(state0, mu))
    if ecc < 1
        # initial guess for circular or elliptical case
        x0 = alpha * sqrt(mu) * (t - t0)
    else
        # initial guess for parabola or hyperbola
        x0 = sqrt(mu) * (t - t0) / (10 * norm(r0))
    end
    x1 = x0    # FIXME - is this the best syntax...?
    # initialize iteratation
    count = 0
    while count < maxiter
        # evaluate function
        Fun, dF, d2F = kepler_der(x0, alpha, t, t0, mu, norm(r0), sigma0)
        if abs(Fun) < tol
            break
        end
        # Laguerre-Conway iteration
        x1 = x0 - 5 * Fun / (dF + dF / abs(dF) * sqrt(abs(16 * dF^2 - 20 * Fun * d2F)))
        count += 1
        x0 = x1
    end

    # ------------------------------------------ #
    # COMPUTE FINAL POSITION
    # convert back to position
    u0, u1, u2, u3 = universal_functions(x1, alpha)
    r_scal = norm(r0) * u0 + sigma0 * u1 + u2
    sigma = sigma0 * u0 + (1 - alpha * norm(r0)) * u1
    #println("r_scal: $r_scal, sigma: $sigma")

    # get lagrange coefficients
    f, g, fdot, gdot = lagrange_coefficients_der(
        mu,
        alpha,
        norm(r0),
        norm(v0),
        sigma0,
        u0,
        u1,
        u2,
        u3,
        r_scal,
        sigma,
    )
    #println("f: $f, g: $g, fdot: $fdot, gdot: $gdot")
    # create map for state   # --- FIXME! can probably speed this up!
    rmap = hcat(f * Matrix(I, 3, 3), g * Matrix(I, 3, 3))
    vmap = hcat(fdot * Matrix(I, 3, 3), gdot * Matrix(I, 3, 3))  #concatenate((fdot*Matrix(I,3,3) , gdot*Matrix(I,3,3)), axis=1)
    fullmap = vcat(rmap, vmap)
    state1 = fullmap * state0

    # ------------------------------------------ #
    # COMPUTE FINAL STM
    # submatrices M's, dyads from pg.377
    M1 = r0 * transpose(r0) / norm(r0)^2
    M2 = r0 * transpose(v0) / (norm(r0) * norm(v0))
    M3 = v0 * transpose(r0) / (norm(r0) * norm(v0))
    M4 = v0 * transpose(v0) / norm(v0)^2
    # coefficients
    c11, c12, c13, c14, c21, c22, c23, c24, c31, c32, c33, c34, c41, c42, c43, c44 =
        der_stm_coefs(mu, alpha, norm(r0), norm(v0), r_scal, u0, u1, u2, u3, sigma0, sigma)
    # construct STM
    Rtilda = f * Matrix(I, 3, 3) + c11 * M1 + c12 * M2 + c13 * M3 + c14 * M4
    R = g * Matrix(I, 3, 3) + c21 * M1 + c22 * M2 + c23 * M3 + c24 * M4
    Vtilda = fdot * Matrix(I, 3, 3) + c31 * M1 + c32 * M2 + c33 * M3 + c34 * M4
    V = gdot * Matrix(I, 3, 3) + c41 * M1 + c42 * M2 + c43 * M3 + c44 * M4
    stm = vcat(hcat(Rtilda, R), hcat(Vtilda, V))
    return state1, stm
end



"""
    keplerder_nostm(mu::Float64, state0::Vector{Float64}, t0::Float64, t::Float64, tol::Float64, maxiter::Int)

Function computes position at some future time t, without computing STM
Formulation from G. Der formulation.

# Arguments
    `mu::Float64`: gravitational parameter
    `state0::Vector{Float64}`: initial state [x,y,z,vx,vy,vz]
    `t0::Float64`: initial time at state0
    `t::Float64`: final time
    `tol::Float64`: tolerance on Laguerre-correction, suggest 1.e-14
    `maxiter::Int`: max allowed iteration allowed for Laguerre-correction, suggest 10

# Returns
    `(array)`: final state
"""
function keplerder_nostm(
    mu::Float64,
    state0::Vector{Float64},
    t0::Float64,
    t::Float64,
    tol::Float64,
    maxiter::Int,
)
    # ------------------------------------------ #
    # SET-UP PROBLEM
    r0, v0 = state0[1:3], state0[4:6]
    sma = get_semiMajorAxis(state0, mu)
    alpha = 1.0 / sma
    sigma0 = dot(r0, v0) / sqrt(mu)

    # ------------------------------------------ #
    # ITERATION WITH LAGUERRE-CONWAY METHOD
    # initial guess based on eccentricity
    ecc = norm(get_eccentricity(state0, mu))
    if ecc < 1
        # initial guess for circular or elliptical case
        x0 = alpha * sqrt(mu) * (t - t0)
    else
        # initial guess for parabola or hyperbola
        x0 = sqrt(mu) * (t - t0) / (10 * norm(r0))
    end
    # initialize final guess
    x1 = x0
    # initialize iteratation
    count = 0
    Fun = 1.0  # FIXME - need storage
    while count < maxiter
        # evaluate function
        Fun, dF, d2F = kepler_der(x0, alpha, t, t0, mu, norm(r0), sigma0)
        if abs(Fun) < tol
            break
        end
        # Laguerre-Conway iteration
        x1 = x0 - 5 * Fun / (dF + dF / abs(dF) * sqrt(abs(16 * dF^2 - 20 * Fun * d2F)))
        # else
        #     # Newton iteration
        #     x1 = x0 - Fun / dF
        # end
        count += 1
        x0 = x1
    end

    # ------------------------------------------ #
    # COMPUTE FINAL POSITION
    # convert back to position
    u0, u1, u2, u3 = universal_functions(x1, alpha)
    r_scal = norm(r0) * u0 + sigma0 * u1 + u2
    sigma = sigma0 * u0 + (1 - alpha * norm(r0)) * u1

    # get lagrange coefficients
    f, g, fdot, gdot = lagrange_coefficients_der(
        mu,
        alpha,
        norm(r0),
        norm(v0),
        sigma0,
        u0,
        u1,
        u2,
        u3,
        r_scal,
        sigma,
    )
    # create map for state   # --- FIXME! can probably speed this up!
    rmap = hcat(f * Matrix(I, 3, 3), g * Matrix(I, 3, 3))
    vmap = hcat(fdot * Matrix(I, 3, 3), gdot * Matrix(I, 3, 3))
    fullmap = vcat(rmap, vmap)
    state1 = fullmap * state0
    return state1
end
