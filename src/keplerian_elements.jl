"""
Function to get Keplerian Elements
Yuri Shimane, 2021/03/23
"""


"""
funtion to get orbital elements
"""
function get_inclination(state)
    """Function computes inclination in radians from a two-body state vector, in inertially frozen frame

    Args:
        state (array): array of cartesian state in inertially frozen frame

    Returns:
       ::Float64: inclination in radians
    """
    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = cross(r, v)
    # inclination
    inc = acos(h[3] / norm(h))
    return inc
end


function get_raan(state)
    """Function computes RAAN in radians from a two-body state vector, in inertially frozen frame

    Args:
        state (array): array of cartesian state in inertially frozen frame

    Returns:
       ::Float64: RAAN in radians
    """

    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = cross(r, v)
    # normal direction of xy-plane
    zdir = [0, 0, 1]
    ndir = cross(zdir, h)
    # compute RAAN
    Ω = atan(ndir[2], ndir[1])
    if Ω < 0.0
        return 2π + Ω
    else
        return Ω
    end
end


function get_eccentricity(state, mu::Float64)
    """Function computes eccentricity vector from a two-body state vector, in inertially frozen frame

    Args:
        state (array): array of cartesian state in inertially frozen frame
        mu::Float64: two-body mass parameter

    Returns:
        (arr): eccentricity vector
    """

    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = cross(r, v)
    ecc = (1 / mu) * cross(v, h) - r / norm(r)
    return ecc
end


function get_aop(state, mu::Float64)
    """Function computes argument of periapsis in radians from a two-body state vector, in inertially frozen frame.
    If eccentricity is 0, omega = RAAN

    Args:
        state (array): array of cartesian state in inertially frozen frame
        mu::Float64: two-body mass parameter

    Returns:
       ::Float64: argument of periapsis in radians
    """

    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = cross(r, v)
    # normal direction of xy-plane
    zdir = [0, 0, 1]
    ndir = cross(zdir, h)
    # compute eccentricity vector
    ecc = (1 / mu) * cross(v, h) - r / norm(r)
    if norm(ecc) != 0
        # compute argument of periapsis
        if dot(ndir, ecc) / (norm(ndir) * norm(ecc)) == -1.0000000000000002  # repeated error case
            ω = 0.0
        elseif dot(ndir, ecc) / (norm(ndir) * norm(ecc)) <= 1.0
            ω = acos(dot(ndir, ecc) / (norm(ndir) * norm(ecc)))
        else
            ω = 0.0   # stability for 1.0000000000000002: case
        end
        if ecc[3] < 0
            ω = 2π - ω
        end
    else
        ω = 0.0  #get_raan(state)
    end
    return ω
end


function get_trueanomaly(state, mu::Float64)
    """Function computes argument of periapsis in radians from a two-body state vector, in inertially frozen frame

    Args:
        state (array): array of cartesian state in inertially frozen frame
        mu::Float64: two-body mass parameter

    Returns:
       ::Float64: true anomaly in radians
    """
    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = norm(cross(r, v))
    # eccentricity vector
    ecc = get_eccentricity(state, mu)
    # radial velocity
    vr = dot(v, r) / norm(r)
    if dot(r, ecc) / (norm(r) * norm(ecc)) <= 1.0
        θ = atan(h * vr, h^2 / norm(r) - mu)
        #acos(dot(r, ecc) / (norm(r)*norm(ecc)))
    else
        θ = acos(1.0)
    end
    # wrap around 2π
    if θ < 0
        return 2π + θ
    else
        return θ
    end
end


function get_semiMajorAxis(state, mu::Float64)
    """Function computes semi major axis of keplrian orbit

    Args:
        state (array): array of cartesian state in inertially frozen frame
        mu::Float64: two-body mass parameter

    Returns:
       ::Float64: semi-major axis
    """
    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = norm(cross(r, v))
    # eccentricity
    e = norm(get_eccentricity(state, mu))
    # semi-major axis
    if abs(1 - e^2) < 1.e-16
        a = Inf
    else
        a = h^2 / (mu * (1 - e^2))
    end
    return a
end


function get_period(state, mu::Float64)
    """Function computes period of keplerian orbit

    Args:
        state (array): array of cartesian state in inertially frozen frame
        mu::Float64: two-body mass parameter

    Returns:
       ::Float64: period
    """
    a = get_semiMajorAxis(state, mu)
    if a > 0
        period = 2 * pi * sqrt(a^3 / mu)
    else
        period = Inf
    end
    return period
end


function get_fpa(state, mu::Float64)
    """Get flight-path angle, between perpendicular to position vector and velocity vector

    Args:
        `state (array)`: state-vector
        `mu::Float64`: gravitational parameter

    Returns:
       ::Float64: FPA, in radians

    """
    h = norm(cross(state[1:3], state[4:6]))
    vperp = h / norm(state[1:3])
    vr = mu / h * norm(get_eccentricity(state, mu)) * sin(get_trueanomaly(state, mu))
    gamma = atan(vr / vperp)
    return gamma
end


"""
Conversion from Cartesian to Keplerian
"""
function cartesian_to_keplerian(state::Vector, μ::Float64, fictious_vz::Float64 = 1.e-15)
    """Get Keplerian elements from Cartesian state-vector and gravitational parameter
    Returned tuple: ecc, sma, inc, raan, aop, ta

    Args:
        `state::Vector{FLoat64}`: state-vector
        `mu::Float64`: gravitational parameter
        `fictious_vz::Float64`: fictious vz to be injected if state is purely planar

    Returns:
        (tuple): eccentricity, semi-major axis, inclination, raan, argument of periapsis, true anomaly
    """
    if state[3] == 0.0 && state[6] == 0.0
        state[6] = fictious_vz
    end
    ecc = get_eccentricity(state, μ)
    sma = get_semiMajorAxis(state, μ)
    inc = get_inclination(state)
    Ω = get_raan(state)  #mod(get_raan(state), 2π)
    ω = get_aop(state, μ)  #mod(get_aop(state, μ), 2π)
    θ = get_trueanomaly(state, μ)  # mod(get_trueanomaly(state, μ), 2π)
    return [norm(ecc), sma, inc, Ω, ω, θ]
end


"""
Conversion from Keplerian to Cartesian
"""
function keplerian_to_cartesian(
    ecc::Float64,
    sma::Float64,
    inc::Float64,
    Ω::Float64,
    ω::Float64,
    θ::Float64,
    μ::Float64,
)
    # angular momentum
    if ecc <= 1.0
        h = sqrt(sma * μ * (1 - ecc^2))
    else
        h = sqrt(abs(sma) * μ * (ecc^2 - 1))   # hyperbolic case
    end
    # perifocal vector
    x = (h^2 / μ) * (1 / (1 + ecc * cos(θ))) * cos(θ)
    y = (h^2 / μ) * (1 / (1 + ecc * cos(θ))) * sin(θ)
    vx = (μ / h) * (-sin(θ))
    vy = (μ / h) * (ecc + cos(θ))
    rpf = [x, y, 0.0]
    vpf = [vx, vy, 0.0]

    r_rot3 = _perifocal2geocentric(rpf, ω, inc, Ω)
    v_rot3 = _perifocal2geocentric(vpf, ω, inc, Ω)

    # save inertial state
    state_inr = vcat(r_rot3, v_rot3)[:] # cat(r_rot3, v_rot3, dims=(1,1))
    return state_inr
end


"""
    keplerian_to_modifiedEquinoctial(keplerian_elts::Vector, use_ta::Bool=true)

Convert Keplerian elements to Modified Equinoctial Elements
Input in order `[ecc, sma, inc, Ω, ω, θ]`, returned in order `[p, f, g, h, k, L]`
Ref: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
"""
function keplerian_to_modifiedEquinoctial(keplerian_elts::Vector, use_ta::Bool = true)
    # unpack
    ecc, sma, inc, Ω, ω, θ = keplerian_elts
    # convert
    p = sma * (1 - ecc^2)
    f = ecc * cos(ω + Ω)
    g = ecc * sin(ω + Ω)
    h = tan(inc / 2) * cos(Ω)
    k = tan(inc / 2) * sin(Ω)
    if use_ta == true
        L = Ω + ω + θ
    else
        # compute eccentric or hyperbolic anomaly
        if ecc <= 1.0
            EA = 2 * atan(sqrt(1 - ecc) * tan(θ / 2), sqrt(1 + ecc))
        else
            EA = 2 * atan(sqrt(ecc - 1) * tan(θ / 2), sqrt(ecc + 1))
        end
        L = Ω + ω + EA
    end
    return [p, f, g, h, k, L]
end


"""
    modifiedEquinoctial_to_keplerian(mee_elts::Vector)

Convert Keplerian elements to Modified Equinoctial Elements
Input in order `[p, f, g, h, k, L]`, returned in order `[ecc, sma, inc, Ω, ω, θ]`
Ref: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
"""
function modifiedEquinoctial_to_keplerian(mee_elts::Vector)
    # unpack
    p, f, g, h, k, L = mee_elts
    # convert
    sma = p / (1 - f^2 - g^2)
    ecc = sqrt(f^2 + g^2)
    inc = atan(2 * sqrt(h^2 + k^2), 1 - h^2 - k^2)
    ω = atan(g * h - f * k, f * h + g * k)
    Ω = atan(k, h)
    θ = L - atan(g, h)
    return [ecc, sma, inc, Ω, ω, θ]
end
