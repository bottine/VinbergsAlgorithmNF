module VinbergsAlgorithmNF

    using Hecke
    using ToggleableAsserts
    using CoxeterDiagrams
    using Convex, Tulip

    import MathOptInterface
    import AbstractAlgebra
    import LinearAlgebra

    
    # Copied from the ToggleableAsserts source
    module Options

        const options_lock = ReentrantLock()

        @enum LastCoordinateMode IsSquare SameOld
        last_coordinate_mode() = IsSquare
        function set_last_coordinate_mode(v::LastCoordinateMode)
            lock(options_lock) do
                @eval Options last_coordinate_mode() = $v
            end
        end

        @enum ConjugatesBoundCheckMode Exact FloatyPlusExactCheck OnlyFloaty
        conjugates_bound_check_mode() = OnlyFloaty 
        function set_conjugates_bound_check_mode(v::ConjugatesBoundCheckMode)
            lock(options_lock) do
                @eval Options conjugates_bound_check_mode() = $v
            end
        end

        lp_precision() = 1024
        function set_lp_precision(v::Int)
            lock(options_lock) do
                @eval Options lp_precision() = $v
            end
        end

        export set_last_coordinate_mode, last_coordinate_mode,
               set_conjugates_bound_check_mode, conjugates_bound_check_mode,
               set_lp_precision,lp_precision
    end



    include("util.jl")
    include("fields.jl")
    include("quad_forms.jl")
    include("vinbergs_algo.jl")
    include("some_lattices.jl")

    export VinbergData,roots_at_distance_zero,cone_roots

end # module
