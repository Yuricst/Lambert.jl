"""
Test using Lambert solver
"""

using Plots
using DifferentialEquations

gr()

psuh!(LOAD_PATH, "../")
using lambert


# initial and final condition
r1vec = [0.79, 0.0, 0.0]
r2vec = [-0.6, -0.17, 0.005]
mu =  1.0
m = 0
tofs = LinRange(1.2, 10.0, 20)

# prepare plot
ptraj = plot(
	set_aspect=:equal, 
	legend=false, 
	frame_style=:box, 
	xlim=[-1.0, 1.0]
)
scatter!(ptraj, [r1vec[1]], [r1vec[2]], color=:green)
scatter!(ptraj, [r2vec[1]], [r2vec[2]], color=:red)

for tof in tofs

	res = lambert_fast(r1vec, r2vec, tof, m, mu)
	pretty(res)

	# plot if exitflag is 1
	if res.exitflag == 1
		# construct state-vector
		sv0 = vcat(r1vec, res.v1)[:]
		svf = vcat(r2vec, res.v2)[:]

		# construct problem
		trajprob_fwd = TwoBodyProblem(0.0, tof, sv0, mu)
		trajprob_bck = TwoBodyProblem(0.0, -tof, svf, mu)

		# solve
		sol_fwd = propagate(trajprob_fwd, Tsit5())
		sol_bck = propagate(trajprob_bck, Tsit5())

		# append plot
		plot!(ptraj, sol_fwd, vars=(1,2))
	end
end

display(ptraj)
println("Done!")


