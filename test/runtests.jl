"""Test sets"""

push!(LOAD_PATH, joinpath(@__DIR__, "../src/"))

using Lambert
using Test

@testset "lambert" begin
    # initial and final condition
    r1vec = [0.79, 0.0, 0.0]
    r2vec = [-0.6, -0.17, 0.015]
    mu =  1.0
    m = 0
    tofs = LinRange(1.2, 10.0, 20)

    for tof in tofs

        res = Lambert.lambert_fast(r1vec, r2vec, tof, m, mu)
		teval = LinRange(0.0, tof, 2)
		traj = Lambert.keplerder_nostm(mu, [r1vec; res.v1], 0.0, teval)

        @test res.exitflag == 1
        @test res.r1 ≈ r1vec atol = 1e-11
        @test res.r2 ≈ r2vec atol = 1e-11
        @test traj[1:3,end] ≈ r2vec atol = 1e-11
    end
end