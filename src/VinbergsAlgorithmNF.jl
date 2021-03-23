module VinbergsAlgorithmNF

    using Hecke
    using ToggleableAsserts
    using CoxeterDiagrams
    using Convex, COSMO, Cbc, Tulip

    import MathOptInterface
    import AbstractAlgebra
    import LinearAlgebra

    include("util.jl")
    include("fields.jl")
    include("quad_forms.jl")
    include("vinbergs_algo.jl")
    include("some_lattices.jl")

    export VinbergData,roots_at_distance_zero,cone_roots

end # module
