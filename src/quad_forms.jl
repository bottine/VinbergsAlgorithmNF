"""
    is_feasible(field,matrix)

    Check that the quadratic form defined by the matrix `matrix` has signature ``(n,1)`` and that all its (non-trivial) Galois conjugates have signature ``(n+1,0)`` (Equivalently, all are definite, the original is hyperbolic and the (non-trivial) conjugates are positive definite).
"""
function is_feasible(V::Hecke.QuadSpace)
    
    P = infinite_places(V.K)
    n = dim(V)

    return signature(V,P[1]) == (dim(V)-1,1,0) && all(signature(V,p) == (dim(V),0,0) for p in P[2:end])
end


function signature(
    V::Hecke.QuadSpace,
    P::InfPlc,
)
    diag = diagonal(V) 
    filter!(≠(0),diag)
    sig =  (count([ispositive(d,P) for d in diag]),count([isnegative(d,P) for d in diag]),dim(V)-length(diag))

    @assert sum(sig) == dim(V)
    return sig
end



function diag_signature(field,diagonal,place)
    non_zeros = filter(≠(0),diagonal)
    (count([ispositive(field(d),place) for d in non_zeros]),count([isnegative(field(d),place) for d in non_zeros]),length(diagonal)-length(non_zeros))
end

function is_diago_and_feasible(field,matrix)
    
    (m,n) = size(matrix)
    
    m ≠ n && return false
    !isdiagonal(matrix) && return false

    diagonal = [matrix[i,i] for i in 1:n]
   
    P = infinite_places(field)

    diag_signature(field,diagonal,P[1]) ≠ (n-1,1,0) && return false
    any(diag_signature(field,diagonal,p) ≠ (n,0,0) for p in P[2:end]) && return false

    # the diagonal is ordered increasingly (this is not really necessary)
    for i in 1:n-1
        diagonal[i] > diagonal[i+1] && return false 
    end

    return true

end





function is_integral(space,ring,vector)
    field = space.K
    
    for c in vector
        if !Hecke.in(field(c),ring)
            return false
        end
    end
    return true
end

function has_positive_norm_squared(space,ring,vector,norm_sq)
    
    @toggled_assert norm_sq == norm_squared(space,vector) "Precomputed norm² must equal actual norm²."
    
    if norm_sq ≤ 0 
        return false
    end
    return true
end

function has_positive_norm_squared(space,ring,vector)
    norm_sq = Hecke.inner_product(space,vector,vector)
    return has_positive_norm_squared(space,ring,vector,norm_sq)
end

function is_primitive(space,ring,vector)
    return isunit(ideal_gcd(ring,vector))
end

function crystallographic_condition(space,ring,vector,norm_sq)
    for b in eachcol(LinearAlgebra.I(length(vector)))
        if  !divides(norm_sq,2*Hecke.inner_product(space,collect(b),vector),ring)
            return false
        end
    end
    return true
end

function crystallographic_condition(space,ring,vector)
    norm_sq = Hecke.inner_product(space,vector,vector)
    return crystallographic_condition(space,ring,vector,norm_sq)
end

function is_root(
    space::Hecke.QuadSpace,
    ring::NfAbsOrd,
    vector::Vector,
    norm_sq,
)

    @toggled_assert norm_sq == norm_squared(space,vector) "Precomputed norm² must equal actual norm²."

    @debug "is_root($vector)"

    !has_positive_norm_squared(space,ring,vector,norm_sq) && return false
    
    @debug "✓ positive length"

    !is_integral(space,ring,vector) && return false
    
    @debug "✓ integral"

    !is_primitive(space,ring,vector) && return false
    
    @debug "✓ primitive"

    !crystallographic_condition(space,ring,vector,norm_sq) && return false 
    
    @debug "✓ crystallographic"

    #Not gonna work since it is only up to squared units: neeed to use the morphisms constructed in possible_root_norms_up_to_squared_units
    #@assert Hecke.inner_product(space,vector,vector) ∈ possible_root_norms_up_to_squared_units(space,ring)

    return true

end

function is_root(
    space::Hecke.QuadSpace,
    ring::NfAbsOrd,
    vector::Vector,
)
    norm_sq = norm_squared(space,vector)
    return is_root(space,ring,vector,norm_sq)
end

function possible_root_norms_squared_up_to_squared_units(
    ring,
    field,
    space,
)

    units, morphism_units = unit_group(ring)
    twice_units, morphism_twice_units = quo(units, 2) # quo(U,2) yields Q = U/U² and the quotient morphism mQ: U -> Q
    representatives_of_units_up_to_squares = [ morphism_units(preimage(morphism_twice_units, q)) for q in twice_units]

    # is it true in general that the lengths divide twice the last invariant factor?
    # PROOF??? TODO

    # last invariant factor is det(G) / gcd(all elements of the cofactor matrix)
    # We use ideals to compute the gcd
    gram = space.gram 
    cofactors = det(gram) * inv(gram)
    gcd_cofactors = ideal_gcd(ring,collect(cofactors))

    # this is the ideal ⟨2*last_invariant_factor⟩ = ⟨2*det(gram)/gcd(cofactors(gram))⟩
    # constructed by taking the ideal I :=⟨gcd(cofactors(gram))⟩
    # and the ideal                   J := ⟨det(gram)⟩
    # then taking the ideal (I:J) = \{x: xI ⊆ J\}
    # then multiplying by the ideal ⟨2⟩
    # Probably one can do something cleaner
    #
    # TODO: probably can get the last invariant factor from Smith Normal Form (should exist in Hecke)
    twice_last_invariant_factor = ideal(ring,2)*colon(ideal(ring,ring(det(gram))),gcd_cofactors)
    twice_last_invariant_factor_factors = Hecke.factor(twice_last_invariant_factor)
    @assert all(Hecke.isprincipal(idl)[1] for idl in keys(twice_last_invariant_factor_factors))

    unit_factors_of_root_lengths = Dict([u=>1 for u in representatives_of_units_up_to_squares])
    prime_factors_of_root_lengths = Dict([Hecke.isprincipal(idl)[2] => mul for (idl,mul) in twice_last_invariant_factor_factors])
    all_factors_of_root_lengths = merge(prime_factors_of_root_lengths, unit_factors_of_root_lengths)

    all_root_norms = 
    filter(
        l -> istotally_positive(field(l)),
        products(all_factors_of_root_lengths)
    )
    return [ring(l) for l in unique(all_root_norms)]

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

function Coxeter_coeff(space, ring, r₁, r₂)
    
    @toggled_assert is_root(space, ring, r₁) && is_root(space, ring, r₂) "The elements must be roots"

    angle = Hecke.inner_product(space,r₁,r₂)
    cos² = approx(angle^2//(Hecke.inner_product(space,r₁,r₁)*Hecke.inner_product(space,r₂,r₂)))
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

        m = round(Int, π/acos(√cos²))

        return m
    end

end

get_Coxeter_matrix(space, ring, roots) = reduce(hcat,[[Coxeter_coeff(space, ring, r₁,r₂) for r₁ in roots] for r₂ in roots])

