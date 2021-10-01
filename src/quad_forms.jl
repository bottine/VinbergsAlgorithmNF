"""
    is_feasible(field,matrix)

    Check that the quadratic form defined by the matrix `matrix` has signature ``(n,1)`` and that all its (non-trivial) Galois conjugates have signature ``(n+1,0)`` (Equivalently, all are definite, the original is hyperbolic and the (non-trivial) conjugates are positive definite).
"""
function is_feasible(field, mat)
    is_feasible(Hecke.quadratic_space(field, matrix(field, mat)))
end

function is_feasible(V::Hecke.QuadSpace)
    
    P = infinite_places(V.K)
    n = dim(V)

    return signature(V, P[1]) == (dim(V) - 1, 1, 0) && all(signature(V, p) == (dim(V), 0, 0) for p in P[2:end])
end


function signature(
    V::Hecke.QuadSpace,
    P::InfPlc,
)
    diag = diagonal(V) 
    filter!(≠(0), diag)
    sig =  (count([ispositive(d, P) for d in diag]), count([isnegative(d, P) for d in diag]), dim(V) - length(diag))

    @assert sum(sig) == dim(V)
    return sig
end



function diag_signature(field, diagonal, place)
    non_zeros = filter(≠(0), diagonal)
    (count([ispositive(field(d), place) for d in non_zeros]), count([isnegative(field(d), place) for d in non_zeros]), length(diagonal) - length(non_zeros))
end

function is_diago_and_feasible(field, matrix)
    
    (m, n) = size(matrix)
    
    m ≠ n && return false
    !isdiagonal(matrix) && return false

    diagonal = [matrix[i,i] for i in 1:n]
   
    P = infinite_places(field)

    diag_signature(field, diagonal, P[1]) ≠ (n - 1, 1, 0) && return false
    any(diag_signature(field, diagonal, p) ≠ (n, 0, 0) for p in P[2:end]) && return false

    # the diagonal is ordered increasingly (this is not really necessary)
    for i in 1:n - 1
        diagonal[i] > diagonal[i + 1] && return false 
    end

    return true

end




function is_integral(space, ring, vector)
    field = space.K
    
    for c in vector
        if !Hecke.in(field(c), ring)
            return false
        end
    end
    return true
end

function has_positive_norm_squared(space, ring, vector, norm_sq)
    
    @toggled_assert norm_sq == norm_squared(space, vector) "Precomputed norm² must equal actual norm²."
    
    if norm_sq ≤ 0 
        return false
    end
    return true
end

function has_positive_norm_squared(space, ring, vector)
    norm_sq = Hecke.inner_product(space, vector, vector)
    return has_positive_norm_squared(space, ring, vector, norm_sq)
end

function is_primitive(space, ring, vector)
    return isunit(ideal_gcd(ring, vector))
end

function crystallographic_condition(space, ring, vector, norm_sq)
    for b in eachcol(LinearAlgebra.I(length(vector)))
        if !divides(norm_sq, 2 * Hecke.inner_product(space, collect(b), vector), ring)
            return false
        end
    end
    return true
end

function crystallographic_condition(space, ring, vector)
    norm_sq = Hecke.inner_product(space, vector, vector)
    return crystallographic_condition(space, ring, vector, norm_sq)
end

function is_root(
    space::Hecke.QuadSpace,
    ring::NfAbsOrd,
    vector::Vector,
    norm_sq,
)

    @toggled_assert norm_sq == norm_squared(space, vector) "Precomputed norm² must equal actual norm²."

    @debug "is_root($vector)"

    !has_positive_norm_squared(space, ring, vector, norm_sq) && return false
    
    @debug "✓ positive length"

    !is_integral(space, ring, vector) && return false
    
    @debug "✓ integral"

    !is_primitive(space, ring, vector) && return false
    
    @debug "✓ primitive"

    !crystallographic_condition(space, ring, vector, norm_sq) && return false 
    
    @debug "✓ crystallographic"

    # Not gonna work since it is only up to squared units: neeed to use the morphisms constructed in possible_root_norms_up_to_squared_units
    # @assert Hecke.inner_product(space,vector,vector) ∈ possible_root_norms_up_to_squared_units(space,ring)

    return true

end

function is_root(
    space::Hecke.QuadSpace,
    ring::NfAbsOrd,
    vector::Vector,
)
    norm_sq = norm_squared(space, vector)
    return is_root(space, ring, vector, norm_sq)
end

# Overload some of Hecke's function for SNF to work 
# Thanks Tommy Hofmann
Hecke.canonical_unit(x::NfOrdElem) = one(parent(x))
Hecke.gcdx(x::NfOrdElem, y::NfOrdElem) = Hecke._gcdx(x, y)
Hecke.div(x::NfOrdElem, y::NfOrdElem) = x
# End of the overloads

"""
    possible_root_norms_squared_up_to_squared_units(ring,field,space)

Compute the set of possible primitive root lengths, up to squared units.

"""
function possible_root_norms_squared_up_to_squared_units(
    ring,
    field,
    gram,
)
   
    units, morphism_units = unit_group(ring)
    # quo(U,2) yields Q = U/U² and the quotient morphism mQ: U -> Q
    units_squared, morphism_units_squared = quo(units, 2)
    # take preimages of U/U², i.e. representatives 
    representatives_of_units_up_to_squares = [ morphism_units(preimage(morphism_units_squared, q)) for q in units_squared] 

    lif = snf(matrix(ring, gram.entries))[end,end]

    twice_last_invariant_factor = ideal(ring, 2 * lif)
    twice_last_invariant_factor_factors = Hecke.factor(twice_last_invariant_factor)
    @assert all(Hecke.isprincipal(idl)[1] for idl in keys(twice_last_invariant_factor_factors))

    unit_factors_of_root_lengths = Dict([u => 1 for u in representatives_of_units_up_to_squares])
    prime_factors_of_root_lengths = Dict([Hecke.isprincipal(idl)[2] => mul for (idl, mul) in twice_last_invariant_factor_factors])
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
    for (c₁, c₂) in zip(r₁, r₂)
        c₁ == 0 && c₂ == 0 && continue
        c₁ == 0 && c₂ ≠ 0  && return false
        c₁ ≠ 0 && c₂ == 0  && return false
        !isnothing(ratio) && c₁ * ratio ≠ c₂ && return false
        
        isnothing(ratio) && (ratio = c₂ // c₁)
    end

    return true

end

function has_primitive_root_of_length(ring, gram, l) 
    # 2Qr = lv iff 
    # (Q⊕0 - 0⊕l*Id) * (r ⊕ v) = 0 
    n = size(gram)[1]

    M = matrix(ring, ring.([gram diagonal_matrix(ring(l), n)])) 

    kerdim, ker = kernel(M)

    return isunit(ring(elem_gcd(ring, ker[1:n,1:kerdim])))
end



function Gram_coeff(space, r₁, r₂)
    return Hecke.inner_product(space, r₁, r₂)^2 // (Hecke.inner_product(space, r₁, r₁) * Hecke.inner_product(space, r₂, r₂))
end

function Coxeter_coeff(space, ring, r₁, r₂)
    

    @toggled_assert is_root(space, ring, r₁) && is_root(space, ring, r₂) "The elements must be roots"

    
    cos² = approx(Gram_coeff(space, r₁, r₂))
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

        m = round(Int, π / acos(√cos²))

        return m
    end

end

function Coxeter_matrix(space, ring, roots) 
    if isempty(roots)
        reshape(Int[], 0, 0)
    else
        return reduce(hcat, [[Coxeter_coeff(space, ring, r₁, r₂) for r₁ in roots] for r₂ in roots])
    end
end
