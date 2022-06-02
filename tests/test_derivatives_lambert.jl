"""
Test with lambert derivatives
"""

using LinearAlgebra

push!(LOAD_PATH, "../")
using Lambert

# settings
r1 = [1,0,0.03]
r2 = [-0.1, 1.2, -0.01]
tof = 2.1
μ = 1.0
cw = false
m = 0

# solve Lambert problem with Jacobians
res, dv_dt, dv_dr = lambert_jacobians(
    r1,
    r2,
    tof,
    m,
    μ,
    cw,
);

println("dv_dt: \n$dv_dt")
println("dv_dr: \n$dv_dr")