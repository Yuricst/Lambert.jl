"""
Trajectory propagation using DifferentialEquations
"""

using DifferentialEquations


abstract type AbstractTrajectoryProblem end

struct TwoBodyProblem <: AbstractTrajectoryProblem
	t0#::Union{Int, Float64}
	tf#::Union{Int, Float64}
	x0#::Vector
	p#::Vector
end


function rhs_twobody!(du,u,p,t)
	# radius
	r = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
	# velocities
	du[1] = u[4]
	du[2] = u[5]
	du[3] = u[6]
	# accelerations
	du[4] = -p[1]/r^3 * u[1]
	du[5] = -p[1]/r^3 * u[2]
	du[6] = -p[1]/r^3 * u[3]
end


function propagate(
	trajprob::TwoBodyProblem, 
	method, 
	reltol::Float64=1.e-11, 
	abstol::Float64=1.e-12,
)
	# construct problem
	prob = ODEProblem(
		rhs_twobody!, 
		trajprob.x0, 
		(trajprob.t0, trajprob.tf), 
		trajprob.p
	)
	# solve
	sol = solve(prob, method, reltol=reltol, abstol=abstol)
	return sol
end