"""
Test spice-based rdv problem
"""

using joptimise
using LinearAlgebra
using Plots
using SPICE
gr()

push!(LOAD_PATH, "../")
using Lambert

# modify paths as necessary!
if Sys.iswindows()
    spice_dir = "C:\\Users\\yshimane3\\Documents\\spice"
elseif Sys.islinux()
    spice_dir = "/mnt/c/Users/yurio/Documents/spice"
elseif Sys.isapple()
    spice_dir = "/Users/yuri/Documents/spice"
end
# get spice kernels
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))

# construct problem
et0_str="2024-01-01T00:00:00.00" 
etf_str="2028-01-01T00:00:00.00"
rendezvous!, lx, ux, ng, lg, ug, visits, canonical_params = construct_rdv2imp_spice_problem(
	[3,4],
	et0_str,
	etf_str,
)
x0 = rand(2)

## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 500,   # 1500 ~ 2500
    "print_level" => 5,
    "tol" => 1e-6
)

xopt, fopt, info = minimize(rendezvous!, x0, ng; 
	lx=lx, 
	ux=ux, 
	lg=lg, 
	ug=ug, 
	solver="ipopt",
	options=ip_options
)

println(info)

# view result
prop_traj, prop_r1, prop_r2 = view_rdv2imp_problem(
	xopt, visits, canonical_params.et0, canonical_params.tstar, canonical_params.vstar)

# prepare plot
ptraj = plot(
	size=(500,500),
	set_aspect=:equal, 
	legend=false, 
	frame_style=:box, 
)

# intial and final position
r1vec = prop_r1[:,1]
r2vec = prop_r2[:,1]
scatter!(ptraj, [r1vec[1]], [r1vec[2]], c=:blue)
scatter!(ptraj, [r2vec[1]], [r2vec[2]], c=:crimson)
plot!(ptraj, [0.0, r1vec[1]], [0.0, r1vec[2]], color=:black)
plot!(ptraj, [0.0, r2vec[1]], [0.0, r2vec[2]], color=:black)

# transfer trajectory
plot!(ptraj, prop_traj[1,:], prop_traj[2,:], linestyle=:dash)
#plot!(ptraj, sol_fwd, vars=(1,2))

# inital and final orbit
plot!(ptraj, prop_r1[1,:], prop_r1[2,:], linestyle=:dash, c=:blue)
plot!(ptraj, prop_r2[1,:], prop_r2[2,:], linestyle=:dash, c=:crimson)

display(ptraj)
println("Done!")


