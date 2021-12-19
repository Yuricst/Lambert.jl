"""
Grid-based analysis of two-burn rendez-vous
"""


using Plots
using DifferentialEquations
using joptimise

gr()

include("lambert_fast.jl")
include("twobody.jl")
include("keplerder.jl")


# initial and final condition
x01 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
x02 = [1.5, 0.0, 0.0004, 0.0, sqrt(1/1.5), 0.0]
mu =  1.0
m = 0

function locate_x1(epoch::Float64)
	# get position and velocity vector at epoch
	xf = keplerder_nostm(mu, x01, 0.0, epoch, 1.e-14, 20)
	return xf[1:3], xf[4:6]
end

function locate_x2(epoch::Float64)
	# get position and velocity vector at epoch
	xf = keplerder_nostm(mu, x02, 0.0, epoch, 1.e-14, 20)
	return xf[1:3], xf[4:6]
end

# bounds on epoch