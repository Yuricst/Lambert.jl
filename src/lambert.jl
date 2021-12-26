module Lambert


using LinearAlgebra
using Printf

include("keplerian_elements.jl")
include("keplerder.jl")
include("lambert_fast.jl")

include("optim/cycler.jl")
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

export construct_mga_problem
export construct_mga1dsm_problem
export construct_cycler2_problem, view_cycler2_problem
export construct_lamblt_problem, view_lamblt_problem

end # module
