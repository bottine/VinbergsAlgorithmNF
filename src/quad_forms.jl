


function is_integral(vector)
    
    for c in vector
        if !isinteger(c)
            return false
        end
    end
    return true
end

times(matrix,u,v) = u' * matrix * v

function has_positive_norm_squared(matrix,vector,norm_sq)
    
    @toggled_assert norm_sq == norm_squared(matrix,vector) "Precomputed norm² must equal actual norm²."
    
    if norm_sq ≤ 0 
        return false
    end
    return true
end

function has_positive_norm_squared(matrix,vector)
    norm_sq = times(matrix,vector,vector)
    return has_positive_norm_squared(matrix,vector,norm_sq)
end

function is_primitive(matrix,vector)
    return abs(gcd(ZZ.(vector))) == 1
end

function crystallographic_condition(matrix,vector,norm_sq)
    for b in eachcol(LinearAlgebra.I(length(vector)))
        if  !divides(norm_sq,2*times(matrix,collect(b),vector))
            return false
        end
    end
    return true
end

function crystallographic_condition(matrix,vector)
    norm_sq = times(matrix,vector,vector)
    return crystallographic_condition(matrix,vector,norm_sq)
end

function is_root(
    matrix,
    vector::Vector,
    norm_sq,
)

    @toggled_assert norm_sq == norm_squared(matrix,vector) "Precomputed norm² must equal actual norm²."

    @debug "is_root($vector)"

    !has_positive_norm_squared(matrix,vector,norm_sq) && return false
    
    @debug "✓ positive length"

    !is_integral(vector) && return false
    
    @debug "✓ integral"

    !is_primitive(matrix,vector) && return false
    
    @debug "✓ primitive"

    !crystallographic_condition(matrix,vector,norm_sq) && return false 
    
    @debug "✓ crystallographic"

    #Not gonna work since it is only up to squared units: neeed to use the morphisms constructed in possible_root_norms_up_to_squared_units
    #@assert times(matrix,vector,vector) ∈ possible_root_norms_up_to_squared_units(matrix)

    return true

end

function is_root(
    matrix,
    vector::Vector,
)
    norm_sq = norm_squared(matrix,vector)
    return is_root(matrix,vector,norm_sq)
end

function colinear(
    r₁::Vector,
    r₂::Vector,
)
    
    @toggled_assert length(r₁) == length(r₂)    
    n = length(r₁)
    @toggled_assert n > 0
    
    ratio = nothing
    for (c₁,c₂) in zip(r₁,r₂)
        c₁ == 0 && c₂ == 0 && continue
        c₁ == 0 && c₂ ≠ 0  && return false
        c₁ ≠ 0 && c₂ == 0  && return false
        !isnothing(ratio) && c₁*ratio ≠ c₂ && return false
        
        isnothing(ratio) && (ratio = c₂//c₁)
    end

    return true

end

function Gram_coeff(matrix, r₁, r₂)
    return times(matrix,r₁,r₂)^2//(times(matrix,r₁,r₁)*times(matrix,r₂,r₂))
end

function Coxeter_coeff(matrix, r₁, r₂)
    

    @toggled_assert is_root(matrix, r₁) && is_root(matrix, r₂) "The elements must be roots"

    
    cos² = Gram_coeff(matrix,r₁,r₂)
    if cos² == 0
        return 2
    elseif cos² == 1
        return 0
    elseif cos² > 1
        return 1
    else

        #   cos(π/m)² = r₁⋅r₂ / r₁² r₂² 
        # ⇒ cos(π/m) = √(r₁⋅r₂/r₁r₂)
        # ⇒ m = π/acos(√(r₁⋅r₂/r₁r₂))

        # TODO Is this rounding dangerous?
        #      Can the angle be different from a submultiple of π? if yes, how to deal with it?

        m = round(Int, π/acos(√Float64(cos²)))

        return m
    end

end

function Coxeter_matrix(matrix, roots) 
    if isempty(roots)
        reshape(Int[],0,0)
    else
        return reduce(hcat,[[Coxeter_coeff(matrix, r₁,r₂) for r₁ in roots] for r₂ in roots])
    end
end
