
mutable struct VinbergData

    dim::Int
    field::AnticNumberField
    ring::NfAbsOrd{AnticNumberField,nf_elem}
    gram_matrix::AbstractAlgebra.Generic.MatSpaceElem{nf_elem}
    quad_space::Hecke.QuadSpace{AnticNumberField,AbstractAlgebra.Generic.MatSpaceElem{nf_elem}}

    diagonal_basis::Vector{Vector{NfAbsOrdElem{AnticNumberField,nf_elem}}}
    diagonal_values::Vector{NfAbsOrdElem{AnticNumberField,nf_elem}}
    scaling::Vector{NfAbsOrdElem{AnticNumberField,nf_elem}}
    
    diagonal_change::Array{nf_elem,2}
    diagonal_change_inv::Array{nf_elem,2}


    possible_root_norms_squared_up_to_squared_units::Vector{NfAbsOrdElem{AnticNumberField,nf_elem}}
end


function VinbergData(number_field,gram_matrix)

    (n,m) = size(gram_matrix)
    @assert n == m "The Gram gram_matrix must be square."

    ring_of_integers = maximal_order(number_field)
    quad_space = quadratic_space(number_field, matrix(number_field,gram_matrix))
    #quad_lattice = Hecke.lattice(quad_space)
    
    @assert is_feasible(quad_space) "The quadratic form must be feasible sig (n,1) and all conjugates sig (n+1,0)"
    @assert all(number_field(c) ∈ ring_of_integers for c in gram_matrix) "The Gram matrix must have coefficients in the ring of integers."

    diagonal_basis_vecs,diagonal_values,scaling = diagonalize_and_get_scaling(gram_matrix,ring_of_integers,number_field)
    negative_vector_index = filter(x-> x[2]<0, collect(enumerate(diagonal_values)))[1][1]
    
    if negative_vector_index ≠ 1
        diagonal_basis_vecs[1],diagonal_basis_vecs[negative_vector_index] = diagonal_basis_vecs[negative_vector_index],diagonal_basis_vecs[1]
        diagonal_values[1],diagonal_values[negative_vector_index] = diagonal_values[negative_vector_index],diagonal_values[1]
        scaling[1],scaling[negative_vector_index] = scaling[negative_vector_index],scaling[1]
    end

    negative_vector_index = filter(x-> x[2]<0, collect(enumerate(diagonal_values)))[1][1]
    @assert negative_vector_index == 1
    #@assert is_diago_and_feasible(number_field,gram_matrix) "The Gram matrix must be feasible, diagonal and its diagonal must be increasing."

    return VinbergData(
        n,
        number_field,
        ring_of_integers,
        matrix(number_field,gram_matrix),
        quad_space,
        diagonal_basis_vecs,
        diagonal_values,
        scaling,
        Matrix(matrix(number_field,diagonal_basis_vecs)),
        Matrix(inv(matrix(number_field,diagonal_basis_vecs))),
        ring_of_integers.(possible_root_norms_squared_up_to_squared_units(ring_of_integers, number_field, quad_space)))

end

function VinbergData(gram_matrix)
    K = parent(gram_matrix[1:1])
    VinbergData(K,gram_matrix)
end


diag(vd::VinbergData) = [vd.gram_matrix[i,i] for i in 1:vd.dim]

function to_diag_rep(vd,vec)
    # vec is in canonical coordinattes
    vd.diagonal_change*vec
end

function to_can_rep(vd,vec)
    # vec is in diagonal coordinates
    vd.diagonal_change_inv*vec
end

function basepoint(vd::VinbergData)
    @assert vd.diagonal_values[1] < 0 "sanity check"

    return vd.field.(vd.diagonal_basis[1])
end

times(quad_space::Hecke.QuadSpace,u,v) = Hecke.inner_product(quad_space,u,v)
times(vd::VinbergData,u,v) = times(vd.quad_space,u,v)

norm_squared(quad_space::Hecke.QuadSpace,u) = times(quad_space,u,u)
norm_squared(vd::VinbergData,u) = times(vd.quad_space,u,u)

function fake_dist_to_basepoint(vd,u)
    
    return  times(vd,u,basepoint(vd))^2//norm_squared(vd,u)
end





LeastKByRootNormSquared = Dict{
    NfAbsOrdElem{AnticNumberField,nf_elem},
    Tuple{
        NfAbsOrdElem{AnticNumberField,nf_elem},
        Vector{NfAbsOrdElem{AnticNumberField,nf_elem}}
    }
}



function enumerate_k(vd::VinbergData,l,k_min,k_max)
    
    #@info "enumerate_k (l=$l, k_min=$k_min, k_max=$k_max)"

    ring = vd.ring
    field = vd.field

    α = vd.diagonal_values[1]
    s = vd.scaling[1]

    @assert α < 0 

    M = approx_sum_at_places(field(l*s^2)//field(α),first_place_idx=2)
    P = infinite_places(field)
    k_min_squared_approx = approx(field((k_min*s)^2))
    k_max_squared_approx = approx(field((k_max*s)^2))

    upscaled_candidates = short_t2_elems(ring, k_min_squared_approx-1 , k_max_squared_approx  + M + 1)
    # short_t2_elems only spits out non-zero stuff, and either k or -k, but not both (for any k)

    map!(abs, upscaled_candidates, upscaled_candidates) # Ensures all are positive  
    
    candidates::Vector{nf_elem} = map(k -> k.elem_in_nf//s, upscaled_candidates) # Correct the scaling  

    #Only k>0 # this should not be needed
    # We only need k>0 because k=0 corresponds to hyperplane containing the basepoint, treated elsewhere
    @toggled_assert all(>(0),candidates)
    
   
    # crystallographic_condition
    filter!(
        k -> divides(l,2*k*α,ring),
        candidates,    
    )
    
    # conjugates are of bounded length
    filter!(
        k -> all(≤(α*k^2,l,p) for p in P[2:end]),
        candidates,
    )
      
    #@info "enumerate_k got $candidates"
    return candidates
end


function next_min_k_for_l(vd,current_k,remaining_k,l)

    k_min,k_max = current_k,current_k+10

    while isempty(remaining_k)
        remaining_k = enumerate_k(vd,l,k_min,k_max)
        filter!(>(current_k),remaining_k)
        sort!(remaining_k)
        k_min,k_max = k_max,k_max+10
    end

    new_k = popfirst!(remaining_k)
    return (new_k,remaining_k)
end

function next_min_k_for_l!(vd,dict,l)

    (current_k,remaining_k) = dict[l] 
    (new_k,new_remaining_k) = next_min_k_for_l(vd,current_k,remaining_k,l)
    dict[l] = (new_k,new_remaining_k)

end


function init_least_k_by_root_norm_squared(vd::VinbergData)

    return Dict(
        [
            vd.ring(l) => next_min_k_for_l(vd,0,[],l)
            for l in vd.possible_root_norms_squared_up_to_squared_units
        ]
    )
end


function next_min_pair(vd,dict)
    field = vd.field
    min_pair = nothing

    # that's the value we want to minimize (related to the distance)
    val(k,l) = field(k^2)//field(l)
    

    for (l,(k,remaining)) in dict 
        if isnothing(min_pair) || val(k,l) < val(min_pair...)
            min_pair = (k,l)
        end
    end
    return min_pair
end

function next_min_pair!(vd,dict)
    
    min_pair = next_min_pair(vd,dict)
    next_min_k_for_l!(vd,dict,min_pair[2])
    return min_pair
end

function extend_root_stem(vd::VinbergData,stem,root_length,bounds=[];t2_cache=nothing)
 
    #@info "extend_root_stem $stem $root_length"

    if t2_cache === nothing
        t2_cache = BoundedT2ElemsCache(vd.ring) 
    end

    # helper function; checks that the vector v is all zeros starting from index i
    zero_from(v,i) = all(==(0),v[i:end])
    
    j = length(stem) + 1
    
    # helper function: update the bounds given by letting the next coeff be k
    bounds_updated(k) = [(r,b - vd.diagonal_values[j]*k*r[j]) for (r,b) in bounds]
 
    #@info "stem is $stem"
    #@info "bounds are:"
    #@info "$([(r[j:end],b) for (r,b) in bounds])"
   
     
    if any( bound < 0 && zero_from(root,j) for (root, bound) in bounds)
        return Vector{Vector{nf_elem}}()
    end
    

    field = vd.field
    ring = vd.ring
    space = vd.quad_space
    P = infinite_places(field)


    # stem = [k₀,…,k_j]

    l = root_length
    #@info tab * "extend_root_stem($stem, $root_length)"

    if j == vd.dim + 1
        

        as_lattice_element = sum([stem[i] .* vd.diagonal_basis[i] for i in 1:vd.dim])
        
        #By the case j== vd.dim, the length should already be == l 
        @toggled_assert times(vd,as_lattice_element,as_lattice_element) == l 

        if is_integral(space, ring, as_lattice_element) && is_root(space,ring,field.(as_lattice_element),l) 
            # this is partially redundant since the crystallographic condition should hold already, as the integralness
            #@info tab * "and it's a root of length $l"
            return Vector{Vector{NfAbsOrdElem}}([ring.(as_lattice_element)])
        else
            #@info tab * "and it's bad (length is $(Hecke.inner_product(vd.quad_space,stem,stem)))"
            return Vector{Vector{NfAbsOrdElem}}()
        end

    elseif j == vd.dim
      
        # If k₁,…,k_{j-1} are chosen and k_j is the last one that needs to be found.
        # If (α₁,…,α_j) is the diagonalized inner product, and r = (k₁,…,k_j) the root to be found, of normed² == l, then
        # We have
        # 
        #   ∑_{i=1}^{j-1} k_i^2α_i + k_j^2α_j == normed²(r) == l 
        #
        # Which means k_j^2 = (l - ∑_{i=1}^{j-1}k_i^2α_i)/α_j.

        issquare,square_root = Hecke.issquare((l - sum([stem[i]^2*vd.diagonal_values[i] for i in 1:j-1]))//vd.diagonal_values[j])
        
        if issquare
            return vcat([extend_root_stem(vd,vcat(stem,[k]),root_length,bounds_updated(k)) for k in unique([square_root,-square_root])]...)            
        else
            return Vector{Vector{NfAbsOrdElem}}()
        end

    else
        
        #@info tab * "stem is not complete"
        

        α = vd.diagonal_values[j]
        s = vd.scaling[j]
        S_j = l - sum([vd.diagonal_values[i]*stem[i]^2 for i in 1:length(stem)]) 
        
        #upscaled_candidates_k_j = short_t2_elems(vd.ring, 0,approx_sum_at_places(field(S_j*s^2)//field(α),first_place_idx=1)+1) # TODO the +1 here is because of inexact computations --> can we make everything exact? --> yes in this case since here approx_sum_at_places ranges over all places, so probably just a trace or an exact t2 computation
        
        conjugates_of_bounded_length(k) = all(≤(field(α*(k//s.elem_in_nf)^2),field(S_j),p) for p in P)
        crystal(k) = divides(l,2*(k//s.elem_in_nf)*α,ring)

        upscaled_candidates_k_j = bounded_t2_elems(
            vd.ring, 
            approx_sum_at_places(field(S_j*s^2)//field(α),first_place_idx=1)+1, 
            t2_cache, 
            [conjugates_of_bounded_length,crystal,],
        )

        candidates_k_j::Vector{nf_elem} = map(x -> x.elem_in_nf//s, upscaled_candidates_k_j)
            @assert candidates_k_j[1] == 0
        
        # bounded length for all conjugates
        #=filter!(
            k-> all(≤(field(α*k^2),field(S_j),p) for p in P),
            candidates_k_j
        )=#
       

        # crystallographic condition
        #=filter!(
            k -> divides(l,2*k*α,ring),
            candidates_k_j,    
        )=#
       
        # add the opposites (short_t2_elems only give one of ±k), and both the crystallographic condition and the bound on length don't distinguish +k from -k, 
        # so we check them for either representatives, and add the others afterwrds
        candidates_k_j = vcat(candidates_k_j, .- candidates_k_j[2:end])
        # end-1 since the end is zero, which we don't want twice
        

        return vcat([extend_root_stem(vd,vcat(stem,[k]),root_length,bounds_updated(k),t2_cache=t2_cache) for k in candidates_k_j]...)
    end
    
end


function roots_at_distance_zero(vd::VinbergData)
    stems = [([0],l) for l in vd.possible_root_norms_squared_up_to_squared_units] 
    
    return vcat([extend_root_stem(vd,stem...) for stem in stems]...)
end

function cone_roots(vd,roots_at_distance_zero)

    @warn "Cone roots computation are approximative ⇒ double check the results by hand."
    roots_at_distance_zero = [vd.field.(root) for root in roots_at_distance_zero]

    # We put first the roots with integer coordinates to maximize the chance of having them in the cone roots
    # It's not necessary but easier to analyze the output and compare with rgug then
    integer_roots = filter!(r -> all(isinteger,vd.field.(r)), roots_at_distance_zero) 
    sort!(integer_roots)
    prepend!(integer_roots,roots_at_distance_zero)

    
    cone_roots = Vector{Vector{nf_elem}}()
    @debug "starting with $(length(roots_at_distance_zero)) at dist zero"


    for r in roots_at_distance_zero
        @debug "looking at $r"
        if  all((-1)*r ≠ cr for cr in cone_roots)
            @debug "so far so good"
            if is_necessary_halfspace(cone_roots,-vd.gram_matrix.entries*r)
                @debug "degeneration"
                push!(cone_roots,r)
            end
        
        end
        @debug "have $(length(cone_roots)) cone roots" 
    end

    cone_roots = drop_redundant_halfspaces(cone_roots)
    @info "have $(length(cone_roots)) cone roots"
    for r in cone_roots
        @info r
    end
    return cone_roots

end

function cone_roots(vd)
    cone_roots(vd,roots_at_distance_zero(vd))
end

function roots_for_pair(vd,pair,prev_roots;t2_cache=nothing)
    
    if t2_cache === nothing
        t2_cache = BoundedT2ElemsCache(vd.ring) 
    end

    (k,l) = pair

    prev_roots_as_diagonals = [vd.diagonal_change*prev_root for prev_root in prev_roots]
    #@info "roots_for_pair($pair,$prev_roots)"
    roots = extend_root_stem(vd,[k],l,[(prev_root,-k*vd.diagonal_values[1]*prev_root[1]) for prev_root in prev_roots_as_diagonals],t2_cache=t2_cache)
    
    #@toggled_assert all(times(vd,root,prev) ≤ 0 for root in roots for prev in prev_roots)  "All angles should be good"
    # just in case
    filter!(root -> all(times(vd,root,prev) ≤ 0 for prev in prev_roots),roots)
    @toggled_assert all(is_root(vd.quad_space,vd.ring,root) for root in roots) "All outputs of extend_root_stem must be roots"
    @toggled_assert all(norm_squared(vd,root) == l for root in roots) "All outputs of extend_root_stem must have correct length"
    @toggled_assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
    @toggled_assert all(r₁ == r₂ || times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots) "Two roots at same distance to basepoint, and both compatible with previous ones, should be compatible with each other: why? ($roots)"
   
    
    return roots
    
end

function roots_for_next_pair!(vd,dict,prev_roots;t2_cache=nothing)

    pair = next_min_pair!(vd,dict)
    @info "next pair is $pair"

    # Here we should assert that the roots in prev_roots are actually closer than the next min pair can give us, otherwise we will get inconsistencies

    return roots_for_pair(vd,pair,prev_roots,t2_cache=t2_cache)
end

function next_n_roots!(vd,prev_roots,dict,das;n=10,t2_cache=nothing)

    if t2_cache===nothing
        t2_cache = BoundedT2ElemsCache(vd.ring) 
    end

    roots = prev_roots
    #Coxeter_matrix = get_Coxeter_matrix(vd.quad_space, vd.ring, prev_roots) 
    new_roots = []
    while n > 0 


        new_roots = roots_for_next_pair!(vd,dict,roots;t2_cache=t2_cache)
        n = n - length(new_roots)

        for root in new_roots 
            @info "Got new root $root"
            extend!(das,[Coxeter_coeff(vd.quad_space, vd.ring, old,root) for old in roots])
            push!(roots,vd.field.(root))
        end


        if is_finite_volume(das)
            return (true,(roots,dict,das))
        end
    end

    return (false,(roots,dict,das))
end

function next_n_roots!(
    vd::VinbergData,
    prev_roots::Vector;
    n=10,
    t2_cache=nothing
)

        Coxeter_matrix = get_Coxeter_matrix(vd.quad_space, vd.ring, prev_roots) 
        das = build_diagram_and_subs(Coxeter_matrix,vd.dim-1)
        dict = init_least_k_by_root_norm_squared(vd)
    
        return next_n_roots!(vd,prev_roots,dict,das;n=n,t2_cache=t2_cache)
end

function next_n_roots!(
    vd::VinbergData;
    n=10,
)
    roots = cone_roots(vd)

    return next_n_roots!(vd,roots;n=n,t2_cache=nothing)
end

function all_in_one(
    diagonal::Vector{nf_elem},
    n::Int
)
    @assert n>0 && n≤100 "Need a reasonable number of roots to look for."
    @assert length(diagonal) > 2 "Need to have a long enough diagonal."
    @assert all(parent(v) == parent(diagonal[1]) for v in diagonal) "The elements of the diagonal must lie in a common field."
    
    field = parent(diagonal[1])

    return next_n_roots!(VinbergData(field,diagm(diagonal)),n=n)
end
