# load module
using LinearAlgebra
using Printf

## abstract lambert result
abstract type AbstractLambertOut end

struct FastLambertOut <: AbstractLambertOut
	tof::Float64
	r1::Vector
	r2::Vector
	v1::Vector
	v2::Vector
	exitflag::Int
end


function pretty(LambertOut::AbstractLambertOut)
	@printf("TOF: %1.4e\n", LambertOut.tof)
	println("Departure: ")
	@printf("r1: %1.4e %1.4e %1.4e\n", 
		LambertOut.r1[1],
		LambertOut.r1[2], 
		LambertOut.r1[3]
	)
	@printf("v1: %1.4e %1.4e %1.4e\n", 
		LambertOut.v1[1],
		LambertOut.v1[2], 
		LambertOut.v1[3]
	)
	println("Arrival: ")
	@printf("r2: %1.4e %1.4e %1.4e\n", 
		LambertOut.r2[1],
		LambertOut.r2[2], 
		LambertOut.r2[3]
	)
	@printf("v2: %1.4e %1.4e %1.4e\n", 
		LambertOut.v2[1],
		LambertOut.v2[2], 
		LambertOut.v2[3]
	)
end


function sqrt_vector(vec)
	return [sqrt(el) for el in vec]
end


function sin_vec(vec)
	return [sin(el) for el in vec]
end


function acos_safe(val)
	return acos( max(-1.0, min(1.0, val)) )
end



"""
Fast lambert algorithm via Dario Izzo's method

# Args
	- `r1::Vector`: initial position vector
	- `r2::Vector`: final position vector
	- `tof::Float64`: time of flight
	- `m::Int` number of revolutions
	- `μ::Float64`: gravitational parameter
	- `cw::Bool`: use clockwise-motion, default is `false`
	- `tol::Float64`: tolerance on Newton-Raphson
	- `maxiter::Int`: max iteration on Newton-Raphson

# Returns
	`FastLambertOut`: object containing output of lambert problem
"""
function lambert_fast(
	r1_vec::Vector, 
	r2_vec::Vector, 
	tof::Float64, 
	m::Int, 
	μ::Float64,
	cw::Bool=false,
	tol::Float64=1.e-12,
	maxiter::Int=20,
)

    # initialize
	exitflag = 0

    # re-scale w.r.t. initial position vector
    r1 = norm(r1_vec)
   	vstar =  sqrt(μ/r1)
    tstar = r1/vstar

    # re-scale
    r1_vec = r1_vec/r1
    r2_vec = r2_vec/r1 
    tof = tof/tstar

    # relevant geometry parameters (non dimensional)
    mr2_vec = norm(r2_vec)
    r1crossr2 = cross(r1_vec, r2_vec)
    dθ = acos_safe(
    	dot(r1_vec,r2_vec)/mr2_vec
	)
    #dθ = acos( max(-1, min(1, (dot(r1_vec,r2_vec))/mr2_vec) ) )  # FIXME
    #println("Initial dθ: $dθ")

    # decide whether to use the left or right branch (for multi-revolution
    # problems), and the long or short way
    leftbranch = sign(m)
    longway = sign(tof)
    m = abs(m)
    tof = abs(tof)

    # check for direction
    if r1crossr2[3] >= 0
    	cw_is_short = true
    else
    	cw_is_short = false
    end

    if (cw_is_short == false) && (dθ < π)
    	dθ = 2π - dθ
    	longway = -1.0
    elseif (cw_is_short == true) && (dθ > π)
    	dθ = 2π - dθ
    	longway = -1.0
    end
    #println("Corrected dθ: $dθ")

    # derived quantities
    c     = sqrt(1 + mr2_vec^2 - 2*mr2_vec*cos(dθ)) # non-dimensional chord
    s     = (1 + mr2_vec + c)/2                     # non-dimensional semi-perimeter
    a_min = s/2                                    # minimum energy ellipse semi major axis
    Λ     = sqrt(mr2_vec)*cos(dθ/2)/s              # Λ parameter (from BATTIN's book)
    mcr       = norm(r1crossr2) #sqrt(crossprd*crossprd.')  # magnitues thereof
    nrmunit   = r1crossr2/mcr                        # unit vector thereof

    # --------------------------------------------------------- #
    # initialization for Newton-Raphson
    logt = log(tof)

    # single revolution (1 solution)
    if (m == 0)

        # initial values
        inn1 = -0.5233      # first initial guess
        inn2 = +0.3233      # second initial guess
        x1   = log(1 + inn1)# transformed first initial guess
        x2   = log(1 + inn2)# transformed first second guess

        # multiple revolutions (0, 1 or 2 solutions)
        # the returned soltuion depends on the sign of [m]
    else
        # select initial values
        if (leftbranch < 0)
            inn1 = -0.5234 # first initial guess, left branch
            inn2 = -0.2234 # second initial guess, left branch
        else
            inn1 = +0.7234 # first initial guess, right branch
            inn2 = +0.5234 # second initial guess, right branch
        end
        x1 = tan(inn1*π/2) # transformed first initial guess
        x2 = tan(inn2*π/2) # transformed first second guess
    end

    # since (inn1, inn2) < 0, initial estimate is always ellipse
    xx    = [inn1, inn2]
    aa    = [
    	a_min/(1 - xx[1]^2),
    	a_min/(1 - xx[2]^2),
    ]
    bβ = [
    	longway * 2*asin(sqrt(((s-c)/2) / aa[1] )),
    	longway * 2*asin(sqrt(((s-c)/2) / aa[2] )),
	]
    aα = 2*acos_safe(minimum(xx))
    #println("aα: $aα")
	#println("bβ: $bβ")

    # evaluate the time of flight via Lagrange expression
    y12  = [
    	aa[1]*sqrt(aa[1]) * ((aα- sin(aα)) - (bβ[1]-sin(bβ[1])) + 2π*m),
    	aa[2]*sqrt(aa[2]) * ((aα- sin(aα)) - (bβ[2]-sin(bβ[2])) + 2π*m)
    ]
    # aa.*sqrt_vector(aa) .*((aα - sin(aα)) .- (bβ-sin_vec(bβ)) .+ 2π*m)
    #println("y12: $y12")

    # initial estimates for y
    if m == 0
        y1 = log(y12[1]) - logt
        y2 = log(y12[2]) - logt
    else
        y1 = y12[1] - tof
        y2 = y12[2] - tof
    end
    
    # initialize error
    err = Inf
    # storage for x
    x = 0

    # --------------------------------------------------------- #
    # Newton-Raphson iterations
    for iterations = 1:maxiter
        #println("***** Iteration $iterations *****")
        # update xnew
        xnew = (x1*y2 - y1*x2) / (y2-y1)
        # copy-pasted code (for performance)
        if m == 0
        	x = exp(xnew) - 1
    	else 
    		x = atan(xnew)*2/π
    	end
        a = a_min/(1 - x^2)
        #println("a_min: $a_min, x: $x, a: $a")

        if (x < 1) # ellipse
            β = longway * 2*asin(sqrt((s-c)/2/a))
            α = 2*acos_safe(x)
        else # hyperbola
            α = 2*acosh(x)
            β = longway * 2*asinh(sqrt((s-c)/(-2*a)))
        end

        # evaluate tof via Lagrange expression
        if (a > 0)
            tof_iter = a*sqrt(a)*((α - sin(α)) - (β-sin(β)) + 2π*m)
        else
            tof_iter = -a*sqrt(-a)*((sinh(α) - α) - (sinh(β) - β))
        end

        # update y
        if m ==0
        	#println("tof_iter: $tof_iter")
        	ynew = log(tof_iter) - logt
        else 
        	ynew = tof_iter - tof
        end

        # update last two iterations (preventing bouncing)
        x1 = x2
        x2 = xnew
        y1 = y2
        y2 = ynew

        # update error
        err = abs(x1 - xnew)
        #println("tof: $tof_iter ... err: $err")

        # break if error is within tolerance
        if err <= tol
        	exitflag = 1
        	break
        end
    end

    # if failed, return here
    if exitflag == 0
    	return FastLambertOut(
    		tof*T, 
    		r1_vec*r1, 
    		r2_vec*r1, 
    		[NaN, NaN, NaN],
			[NaN, NaN, NaN],
    		exitflag
		)
    end

    # Solution for the semi-major axis
    a = a_min/(1-x^2)

    # calculate ψ
    if (x < 1) # ellipse
        β = longway * 2asin(sqrt((s-c)/2/a))
        # ensure acos is well-defined
        α = 2*acos_safe(x)
        ψ  = (α-β)/2
        η2 = 2a*sin(ψ)^2/s
        η  = sqrt(η2)
    else       # hyperbola
        β = longway * 2asinh(sqrt((c-s)/2/a))
        α = 2acosh(x)
        ψ  = (α-β)/2
        η2 = -2a*sinh(ψ)^2/s
        η  = sqrt(η2)
    end

    # normalized normal vector
    ih = longway * nrmunit

    # compute unit-vectors in r1, r2 directions
    r1n = r1_vec/norm(r1_vec)
    r2n = r2_vec/mr2_vec

    # cross-products of ih with unit r1, r2
    ihcrossr1 = cross(ih,r1_vec)
    ihcrossr2 = cross(ih,r2n)

    # radial and tangential components of v1
    vr1 = 1/η/sqrt(a_min) * (2*Λ*a_min - Λ - x*η)
    vt1 = sqrt(mr2_vec/a_min/η2 * sin(dθ/2)^2)

    # radial and tangential components of v2
    vt2 = vt1/mr2_vec
    vr2 = (vt1 - vt2)/tan(dθ/2) - vr1

    # obtain velocity vectors
    v1 = (vr1*r1n + vt1*ihcrossr1) * vstar
    v2 = (vr2*r2n + vt2*ihcrossr2) * vstar

    # construct output
    return FastLambertOut(tof*tstar, r1_vec*r1, r2_vec*r1, v1, v2, exitflag)
end

