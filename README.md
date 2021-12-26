# Lambert.jl
Fast Lambert problem solver in Julia and associated optimization problem formulation. 

### Dependencies

The module has dependencies to
- `LinearAlgebra`, `Printf`

The examples also make use of
- [`joptimise`](https://github.com/Yuricst/joptimise), `Plots`

### Basic Lambert solver

```julia
push!(LOAD_PATH, pwd())  # push load path if necessary

using Lambert

# initial and final condition
r1 = [0.79, 0.0, 0.0]      # initial position vector
r2 = [-0.6, -0.17, 0.005]  # final position vector
mu =  1.0                  # gravitational parameter
m = 0                      # number of revolution(s), >=0

# call Lambert routine
res = lambert_fast(r1vec, r2vec, tof, m, mu)
```

The function returns an `AbstractLambertOut` type, which contains: 

```
res.r1        # initial position vector
res.r2        # final position vector
res.v1        # initial velocity vector
res.v2        # final velocity vector
res.exitflag  # success (1) or fail (0)
```

To display the summary, a practical method is to use

```julia
pretty(res)
```

