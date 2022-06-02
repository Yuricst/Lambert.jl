"""
    lambert_jacobians(
        r1_vec::Vector,
        r2_vec::Vector,
        tof::Float64,
        m::Int,
        μ::Float64,
        cw::Bool = false,
        tol::Float64 = 1.e-12,
        maxiter::Int = 20,
    )

Solve Lambert problem and get its sensitivities. 
Fast lambert algorithm via Dario Izzo's method.
Sensitivities from Arora et al, 2015. JGCD. 
*Partial derivatives of the solution to the lambert boundary value problem*.
Adopted from Python implementation by @kdricemt

# Args
	- `r1_vec::Vector`: initial position vector
	- `r2_vec::Vector`: final position vector
	- `tof::Float64`: time of flight
	- `m::Int` number of revolutions
	- `μ::Float64`: gravitational parameter
	- `cw::Bool`: use clockwise-motion, default is `false`
	- `tol::Float64`: tolerance on Newton-Raphson
	- `maxiter::Int`: max iteration on Newton-Raphson

# Returns
	`FastLambertOut`: object containing output of lambert problem
	`Matrix{Float64}: 6x2 matrix of dv/dt
	`Matrix{Float64}: 6x6 matrix of dv/dr
"""
function lambert_jacobians(
	r1_vec::Vector,
    r2_vec::Vector,
    tof::Float64,
    m::Int,
    μ::Float64,
    cw::Bool = false,
    tol::Float64 = 1.e-12,
    maxiter::Int = 20,
)
	# solve Lambert problem to get velocity vectors
	res = lambert_fast(
	    r1_vec,
	    r2_vec,
	    tof,
	    m,
	    μ,
	    cw,
	    tol,
	    maxiter,
	)
	if res.exitflag == 0
		return res, NaN, NaN
	end

	v1_vec = res.v1
	v2_vec = res.v2

	# define common variables
    r1n = norm(r1_vec)
    r2n = norm(r2_vec)
    v1n = norm(v1_vec)
    v2n = norm(v2_vec)
    k = [0, 0, 1]
    smu = sqrt(μ)

    # true anamaly difference
    cos_dnu = dot(r1_vec, r2_vec)/r1n/r2n
    sin_dnu = sign(dot(cross(r1_vec, r2_vec), k)) * sqrt(1 - cos_dnu^2)

    energy = v1n^2/2 - μ/r1n
    a = -μ/2/energy
    h = cross(r1_vec, v1_vec)
    e = sqrt(1 - norm(h)^2/μ/a)
    p = a * (1 - e^2)

    # "ミッション解析と軌道力学の基礎" p99
    f = 1 - r2n/p * (1 - cos_dnu)
    fdot = sqrt(μ/p) * ((1 - cos_dnu)/sin_dnu) * (1/p * (1 - cos_dnu) - 1/r1n - 1/r2n)
    g = r1n * r2n * sin_dnu/ sqrt(μ * p)
    gdot = 1 - r1n/p * (1 - cos_dnu)

    # kepler STM
    alpha = 1/a
    sig1 = dot(r1_vec, v1_vec) / smu
    sig2 = dot(r2_vec, v2_vec) / smu
    chi = alpha * smu * tof + sig2 - sig1

    # U expressions (eq. 38 - 41 from Arora et al)
    U1 = - r2n * r1n * fdot / smu
    U2 = r1n * (1 - f)
    U3 = smu * (tof - g)
    U4 = U1*U3 - 1/2 * (U2^2 - alpha*U3^2)
    U5 = (chi^3/6 - U3)/alpha
    Cbar = 1/smu * (3 * U5 - chi * U4 - smu * tof * U2)

    # r1 = r1.reshape(3,1)
    # r2 = r2.reshape(3,1)
    # v1 = v1.reshape(3,1)
    # v2 = v2.reshape(3,1)
    dv = v2_vec - v1_vec   # FIXME - check if this is ok for Julia
    dr = r2_vec - r1_vec

    # matrices (eq. 46 - 49 from Arora et al)
    A = r2n/μ * dv * transpose(dv) + 1/r1n^3 * (r1n * (1 - f) * r2_vec*transpose(r1_vec) + Cbar*v2_vec*transpose(r1_vec)) + f*I(3)
    B = r1n/μ * (1- f) * (dr*transpose(v1_vec) - dv*transpose(r1_vec)) + Cbar/μ * v2_vec*transpose(v1_vec) + g*I(3)
    C = -1/r1n^2 * dv*transpose(r1_vec) - 1/r2n^2 * r2_vec*transpose(dv) + fdot*(I(3) - 1/r2n^2 * r2_vec*transpose(r2_vec) + 1/(μ * r2n) * (r2_vec*transpose(v2_vec) - v2_vec*transpose(r2_vec))*r2_vec*transpose(dv)) - μ * Cbar/r2n^3/r1n^3 * r2_vec*transpose(r1_vec)
    D = r1n/μ * dv*transpose(dv) + 1/r2n^3 * (r1n * (1-f) * r2_vec*transpose(r1_vec) - Cbar * r2_vec*transpose(v1_vec)) + gdot * I(3)
    # println("A: $A")
    # println("B: $B")
    # println("C: $C")
    # println("D: $D")

    r1dot = v1_vec
    r2dot = v2_vec
    v1dot = -μ * r1_vec/r1n^3
    v2dot = -μ * r2_vec/r2n^3

    BinvA = inv(B)*A  # FIXME - syntax for Julia?
    CDBinvA = C - D*BinvA
    DBinv = transpose( (transpose(B) \ transpose(D)) )

    # println("BinvA: $BinvA")
    # println("CDBinvA: $CDBinvA")
    # println("DBinv: $DBinv")

    # time

    # from paper -> does not mach result
    # dv1_dt1 = - (v1dot + BinvA @ r1dot)  # 3 x 1
    # dv1_dt2 = CDBinvA @ r2dot
    # dv2_dt1 = CDBinvA @ r1dot
    # dv2_dt2 = v2dot - DBinv @ r2dot

    dv1_dt1 = (v1dot + BinvA*r1dot)  # 3 x 1
    dv1_dt2 = -dv1_dt1
    dv2_dt1 = - CDBinvA*r1dot
    dv2_dt2 = -dv2_dt1

    # construct matrix for time-sensitivities
    dv_dt = [dv1_dt1 dv1_dt2; dv2_dt1 dv2_dt2]  # 6 x 2

    # position and velocity
    dv1_dr1 = - transpose(BinvA)   # 3 x 3
    dv1_dr2 = - transpose(CDBinvA)
    dv2_dr1 = CDBinvA
    dv2_dr2 = DBinv

    # construct matrix for r-sensitivities
    dv_dr = [dv1_dr1 dv1_dr2; dv2_dr1 dv2_dr2]  # 6 x 6

    return res, dv_dt, dv_dr
end