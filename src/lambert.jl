module Lambert


using LinearAlgebra
using Printf
using SPICE

include("util_expressions.jl")
include("keplerian_elements.jl")
include("keplerder.jl")
include("lambert_fast.jl")
include("canonical.jl")

include("optim/cycler.jl")
include("optim/rdv2imp_problem.jl")
include("optim/mga_problem.jl")
include("optim/mga1dsm_problem.jl")
include("optim/lamblt_problem.jl")

export FastLambertOut
export pretty
export propagate_arc
export lambert_fast
export get_inclination,
    get_raan,
    get_eccentricity,
    get_aop,
    get_trueanomaly,
    get_semiMajorAxis,
    get_period,
    get_fpa,
    cartesian_to_keplerian,
    keplerian_to_cartesian,
    keplerian_to_modifiedEquinoctial,
    modifiedEquinoctial_to_keplerian
export keplerder, keplerder_nostm
export construct_rdv2imp_problem,
    construct_rdv2imp_problem_tfcon, view_rdv2imp_problem, construct_rdv2imp_spice_problem
export construct_mga_problem
export construct_mga1dsm_problem, view_mga1dsm_problem
export construct_cycler2_problem, view_cycler2_problem
export construct_lamblt_problem, view_lamblt_problem

end # module
