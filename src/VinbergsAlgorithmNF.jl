module VinbergsAlgorithmNF

    using Hecke
    using ToggleableAsserts
    using CoxeterDiagrams
    using Convex, Tulip

    import MathOptInterface
    import AbstractAlgebra
    import LinearAlgebra

    const LP_PRECISION = 128::Int

    include("util.jl")
    include("fields.jl")
    include("quad_forms.jl")
    include("vinbergs_algo.jl")
    include("some_lattices.jl")

    export VinbergData,roots_at_distance_zero,cone_roots

end # module
