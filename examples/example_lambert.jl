"""
Test using Lambert solver
"""

using GLMakie

include(joinpath(@__DIR__, "../src/Lambert.jl"))

# initial and final condition
r1vec = [0.79, 0.0, 0.0]
r2vec = [-0.6, -0.17, 0.015]
mu =  1.0
m = 0
tofs = LinRange(1.2, 10.0, 20)

# prepare plot
fig = Figure(size=(500,500))
ax = Axis3(fig[1,1], xlabel="x", ylabel="y", zlabel="z")
scatter!(ax, [r1vec[1]], [r1vec[2]], [r1vec[3]], color=:green)
scatter!(ax, [r2vec[1]], [r2vec[2]], [r2vec[3]], color=:red)

for tof in tofs

	res = Lambert.lambert_fast(r1vec, r2vec, tof, m, mu)
	display(res)

	if res.exitflag == 1
		# teval = LinRange(0.0, tof, 100)
		# traj = Lambert.keplerder_nostm(mu, [r1vec; res.v1], 0.0, teval)
		traj = Lambert.propagate_arc(res)
		lines!(ax, traj[1,:], traj[2,:], traj[3,:], color=:blue)
	end
end

display(fig)
println("Done!")


