
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

    # Sort the basis vectors with increasing number of zeroes in their canonical coordinates:
    # This should improve stem enumeration on non-diagonal matrices: stems that define non-integral coordinates are spotted earlier.
    perm_zero = sortperm(diagonal_basis_vecs[2:end],by=(x->count(==(0),x))) .|> (x->x+1)
    diagonal_basis_vecs[2:end] = diagonal_basis_vecs[perm_zero]
    diagonal_values[2:end] = diagonal_values[perm_zero]
    scaling[2:end] = scaling[perm_zero]

    negative_vector_index = filter(x-> x[2]<0, collect(enumerate(diagonal_values)))[1][1]
    @assert negative_vector_index == 1
    #@assert is_diago_and_feasible(number_field,gram_matrix) "The Gram matrix must be feasible, diagonal and its diagonal must be increasing."

    vd =  VinbergData(
        n,
        number_field,
        ring_of_integers,
        matrix(number_field,gram_matrix),
        quad_space,
        diagonal_basis_vecs,
        diagonal_values,
        scaling,
        Matrix(matrix(number_field,diagonal_basis_vecs))',    # notice the transpose here and below. Much frustration took place before I found out I needed those!
        Matrix(inv(matrix(number_field,diagonal_basis_vecs)))',
        ring_of_integers.(possible_root_norms_squared_up_to_squared_units(ring_of_integers, number_field, quad_space)))

    @info "Matrix is $gram_matrix"
    @info "Basepoint is $(basepoint(vd))"
    return vd

end

function VinbergData(gram_matrix)
    K = parent(gram_matrix[1:1])
    VinbergData(K,gram_matrix)
end


diag(vd::VinbergData) = [vd.gram_matrix[i,i] for i in 1:vd.dim]

function to_diag_rep(vd,vec)
    # vec is in canonical coordinattes
    vd.diagonal_change_inv*vec
end

function to_can_rep(vd,vec)
    # vec is in diagonal coordinates
    vd.diagonal_change*vec
end

function basepoint(vd::VinbergData)
    @assert vd.diagonal_values[1] < 0 "sanity check"

    return vd.field.(vd.diagonal_basis[1])
end

times(quad_space::Hecke.QuadSpace,u,v) = Hecke.inner_product(quad_space,u,v)
times(vd::VinbergData,u,v) = times(vd.quad_space,u,v)

norm_squared(quad_space::Hecke.QuadSpace,u) = times(quad_space,u,u)
norm_squared(vd::VinbergData,u) = times(vd.quad_space,u,u)

"""
    fake_dist_to_point(vd,point,root)

If `point` defines an element of ``ℍ^n`` and `root` a hyperplane, compute the fake distance between the point and the hyperplane, satisfying
```math
    \\sinh²(dist(point,hyperplane)) = fake_dist(point,hyperplane)
```
"""
function fake_dist_to_point(vd,point,root)
    @toggled_assert is_root(vd.quad_space,vd.ring,root) "Need a root"
    @toggled_assert norm_squared(vd,point) < 0 "Need a point in hyperbolic space."

    fake_dist = -times(vd,root,point)^2//(norm_squared(vd,point)*norm_squared(vd,root))

    return fake_dist 
end

function fake_dist_to_basepoint(vd,u)

    

    fake_dist =  fake_dist_to_point(vd,basepoint(vd),u)
    @toggled_assert fake_dist == -(to_diag_rep(vd,u)[1]^2)*vd.diagonal_values[1]//norm_squared(vd,u) 

    return fake_dist

end





const LeastKByRootNormSquared = Dict{
    NfAbsOrdElem{AnticNumberField,nf_elem},
    Tuple{
        NfAbsOrdElem{AnticNumberField,nf_elem},
        Vector{NfAbsOrdElem{AnticNumberField,nf_elem}}
    }
}



function enumerate_k(vd::VinbergData,l,k_min,k_max)
    
    @info "enumerate_k (l=$l, k_min=$k_min, k_max=$k_max)"

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


mutable struct BoundedT2ElemsCache
    bounds::Vector{Int}
    elems::Vector{Vector{NfAbsOrdElem}}
end

function BoundedT2ElemsCache(ring)
    return (BoundedT2ElemsCache)(Vector{Int}([0]),Vector{Vector{NfAbsOrdElem}}([[ring(0)]]))
end

function bounded_t2_elems(
    ring,
    t2_bound,
    cache, 
    filters = []
)
    
    int_bound = ceil(Int,t2_bound)
    if int_bound > cache.bounds[end]
        with_margin = ceil(Int,int_bound*1.2)
        push!(cache.bounds,with_margin)
        new_elems = filter(
            x -> x∉cache.elems[end],
            short_t2_elems(ring,cache.bounds[end-1],cache.bounds[end]) .|> abs,
        )
        #sort!(new_elems,rev=true)
        push!(cache.elems, new_elems)
    end
    
    elems = Vector{NfAbsOrdElem}()
    i=1
    while cache.bounds[i] ≤ t2_bound
        append!(
            elems,
            filter(
                x -> all(f(x) for f in filters),
                cache.elems[i],
            )
        )
        i = i+1
        i > length(cache.bounds) && break
    end
        
    if i ≤ length(cache.bounds) && ( cache.bounds[i-1]≠t2_bound ) 
        # only those have to be checked exactly for boundedness
        append!(
            elems,
            filter(
                x -> Hecke.t2(x)≤ t2_bound && all(f(x) for f in filters),
                cache.elems[i]
            )
        )
    end
    return elems 

end

const AffineConstraints = Tuple{
    Vector{Vector{nf_elem}}, # Previous roots
    Vector{nf_elem},         # ≤ bound
    Vector{Int},           # last non non-zero coordinate
}     
const no_constraints = (Vector{Vector{nf_elem}}(),Vector{nf_elem}(),Vector{Int}())

function clearly_inconsistent(ac::AffineConstraints,idx,dim)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac 
    return any( bound < 0 && last_non_zero < idx for (bound,last_non_zero) in zip(c_values,c_last_non_zero_coordinates))   
end

function update_constraints(ac::AffineConstraints,idx::Int,val_at_idx::nf_elem)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac
    
    return (c_vectors,[b - val_at_idx*r[idx] for (r,b) in zip(c_vectors,c_values)],c_last_non_zero_coordinates)

end

function interval_for_k_j(ac::AffineConstraints,j::Int)::Tuple{Union{Nothing,nf_elem},Union{Nothing,nf_elem}}
    applicable = [(vec[j],val) for (vec,val,last) in zip(ac...) if last == j]
    pos = [v//c for (c,v) in applicable if c > 0]
    neg = [v//c for (c,v) in applicable if c < 0]
    return (isempty(neg) ? nothing : maximum(neg), isempty(pos) ? nothing : minimum(pos))
end

function in_interval(
    k,
    interval::Tuple{Union{Nothing,nf_elem},Union{Nothing,nf_elem}}
)
    (lb,ub) = interval
    return (lb === nothing || k ≥ lb) && (ub === nothing || k ≤ ub)
end


function extend_root_stem(
    vd::VinbergData,
    stem::Vector{nf_elem},
    stem_can_rep::Vector{nf_elem},
    root_length,
    root_length_minus_stem_norm_squared,
    constraints::AffineConstraints,
    t2_cache::BoundedT2ElemsCache
)
 
    #@info "extend_root_stem $stem for length =  $root_length"

    # helper function; checks that the vector v is all zeros starting from index i
    
    j = length(stem) + 1

    (c_vectors,c_values,c_last_non_zero_coordinates) = constraints
    
 
    #@info "stem is $stem"
    #@info "bounds are:"
    #@info "$([(r[j:end],b) for (r,b) in bounds])"
   
     
    if clearly_inconsistent(constraints,j,vd.dim) 
        return Vector{Vector{nf_elem}}()
    end
    

    field = vd.field
    ring = vd.ring
    space = vd.quad_space
    P = infinite_places(field)
    l = root_length



    # stem = [k₀,…,k_j]

    #@info tab * "extend_root_stem($stem, $root_length)"

    if j == vd.dim + 1
        
        @toggled_assert stem_can_rep == to_can_rep(vd,stem) "Sanity check. `stem_can_rep` must always be equal to `to_can_rep(vd,stem)`!"
        @toggled_assert times(vd,stem_can_rep,stem_can_rep) == l "Sanity check. If we are here, by the previous case (j==vd.dim) necessarily we have the right length."
        @toggled_assert all(bound ≥ 0 for bound in c_values) "All affine constraints given by previous roots must be satisfied."

        if is_integral(space, ring, stem_can_rep) && is_root(space,ring,field.(stem_can_rep),l) 
            # this is partially redundant since the integralness should already hold?
            #@info tab * "and it's a root of length $l"
            return Vector{Vector{NfAbsOrdElem}}([ring.(stem_can_rep)])
        else
            #@info "isnoroot"
            #@info tab * "and it's bad (length is $(Hecke.inner_product(vd.quad_space,stem,stem)))"
            return Vector{Vector{NfAbsOrdElem}}()
        end

    end

    α_j = vd.diagonal_values[j]
    v_j = vd.diagonal_basis[j]    
    s_j = vd.scaling[j]
    l_j = root_length_minus_stem_norm_squared # l - sum([vd.diagonal_values[i]*stem[i]^2 for i in 1:length(stem)]) # S_j should be computed incrementally since S_{j+1} = S_j - stem[j+1]^2*diag_value[j+1] 

    if j == vd.dim
      
        # If k₁,…,k_{j-1} are chosen and k_j is the last one that needs to be found.
        # If (α₁,…,α_j) is the diagonalized inner product, and r = (k₁,…,k_j) the root to be found, of normed² == l, then
        # We have
        # 
        #   ∑_{i=1}^{j-1} k_i^2α_i + k_j^2α_j == normed²(r) == l 
        #
        # Which means k_j^2 = (l - ∑_{i=1}^{j-1}k_i^2α_i)/α_j.
        issquare,square_root = Hecke.issquare(l_j//α_j)
        
        
        if issquare && divides(l,2*square_root*α_j,ring) # crystal
            return vcat(
                [   
                 extend_root_stem(vd,vcat(stem,[k]),stem_can_rep + k .* v_j,l,l_j - k^2*α_j,update_constraints(constraints,j,k*α_j),t2_cache) 
                    for k in unique([square_root,-square_root])
                ]...
            )            
        else
            return Vector{Vector{NfAbsOrdElem}}()
        end

    else
        
        #@info "$stem is not complete"
        
        # does not work yet
        #interval_k_j = interval_for_k_j(constraints,j)
        #@info "interval $interval_k_j"
        
        
        #
        #       k^2α ≤ l_j at all places
        #    ⇔  (sk)^2 α/s^2 ≤ l_j at all places
        #
        function conjugates_of_bounded_length(sk)
            sat = all(≤(field(α_j*(sk//s_j.elem_in_nf)^2),field(l_j),p) for p in P)
            return sat
        end
        
        #
        # l | 2kα  ⇔  2kα/l ∈ ring  ⇔  sk (2α//ls) ∈ ring 
        #
        function crystal(sk) 
            sat = divides(l,2*(sk//s_j.elem_in_nf)*α_j,ring)
            return sat
        end

        upscaled_candidates_k_j = bounded_t2_elems(
            vd.ring, 
            approx_sum_at_places(field(l_j*s_j^2)//field(α_j),first_place_idx=1)+1, 
            t2_cache, 
            [conjugates_of_bounded_length,crystal,],
        )# TODO the +1 here is because of inexact computations --> can we make everything exact? --> yes in this case since here approx_sum_at_places ranges over all places, so probably just a trace or an exact t2 computation


        candidates_k_j::Vector{nf_elem} = map(x -> x.elem_in_nf//s_j, upscaled_candidates_k_j)
        @assert isempty(candidates_k_j) || candidates_k_j[1] == 0
        
        
        ###################################### The following only useful for non diagonal forms
        ###################################### This should probably be moved inside the bounded_t2_elems call for consistency with the other filters 
        all_zero_on_coord_after(vd,vector_idx,coord_idx) = all(vd.diagonal_basis[l][coord_idx] == 0 for l in vector_idx+1:vd.dim) # this can be precomputed and stored in vd
        integral(k) = all(
            all_zero_on_coord_after(vd,j,idx) ⇒ (stem_can_rep[idx] + k*(v_j[idx]) ∈ ring) 
            for idx in 1:vd.dim
        )
        ######################################
        
        candidates_k_j = vcat(filter(k -> #=in_interval(k,interval_k_j) && =# integral(k),candidates_k_j), filter(k -> #=in_interval(k,interval_k_j) && =# integral(k), .- candidates_k_j[2:end])) # [2:end] since index 1 contains zero, which we don't want twice
        #=for k in candidates_k_j
            if !in_interval(k,interval_k_j)
                @info "$k out of bounds ? $interval_k_j"
                @info "stem: $stem"
                @info "constraints $constraints"
            end
        end=#

        return vcat([extend_root_stem(vd,vcat(stem,[k]),stem_can_rep + k .* v_j,l,l_j - k^2*α_j,update_constraints(constraints,j,k*α_j),t2_cache) for k in candidates_k_j]...)
    end
    
end


function roots_at_distance_zero(vd::VinbergData)

    t2_cache = BoundedT2ElemsCache(vd.ring) 
    zero_stem = nf_elem[vd.field(0)]
    zero_stem_can_rep = nf_elem[vd.field(0) for i in 1:vd.dim]
    
    return vcat([extend_root_stem(vd,zero_stem,zero_stem_can_rep,l,l,no_constraints,t2_cache) for l in vd.possible_root_norms_squared_up_to_squared_units]...)
end

function cone_roots(vd,roots_at_distance_zero)

    @warn "Cone roots computation are approximative ⇒ double check the results by hand."
    roots_at_distance_zero = [vd.field.(root) for root in roots_at_distance_zero]
    @info "starting with $(length(roots_at_distance_zero)) roots at distance zero"
    
    len = length(roots_at_distance_zero)

    # We put first the roots with integer coordinates to maximize the chance of having them in the cone roots
    # It's not necessary but easier to analyze the output and compare with rgug then
    integer_roots = filter(r -> all(isinteger,vd.field.(r)), roots_at_distance_zero)
    non_integer_roots = filter(r -> !all(isinteger,vd.field.(r)), roots_at_distance_zero)
    sort!(integer_roots)
    roots_at_distance_zero = vcat(integer_roots,non_integer_roots)
    @assert len == length(roots_at_distance_zero) "did we drop a root while reorgarizing them?"
    
    cone_roots = Vector{Vector{nf_elem}}()
    @debug "starting with $(length(roots_at_distance_zero)) at dist zero"


    for r in roots_at_distance_zero
        @debug "looking at $r"
        if  all((-1)*r ≠ cr for cr in cone_roots)
            @debug "so far so good"
            if is_necessary_halfspace(vd.gram_matrix.entries,cone_roots,r)
                @debug "degeneration"
                push!(cone_roots,r)
            end
        
        end
        @debug "have $(length(cone_roots)) cone roots" 
    end
    
    #@info "before dropping have $(length(cone_roots)) roots"
    cone_roots = drop_redundant_halfspaces(vd.gram_matrix.entries,cone_roots)
    #@info "have $(length(cone_roots)) cone roots"
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
    fake_dist = -(k^2)*vd.diagonal_values[1].elem_in_nf//(l.elem_in_nf)
   
    ######################## The following only useful for non-diags
    all_zero_on_coord_after(vd,vector_idx,coord_idx) = all(vd.diagonal_basis[l][coord_idx] == 0 for l in vector_idx+1:vd.dim)
    if !all(all_zero_on_coord_after(vd,1,idx) ⇒ (k*(vd.diagonal_basis[1][idx]) ∈ vd.ring) for idx in 1:vd.dim)
        @info "$k would define non-integral coordinates ⇒ dropping it"
        return [] 
    end
    ######################
    
    # Construct the constraints defined by the previous roots
    prev_roots_diag = [to_diag_rep(vd,prev_root) for prev_root in prev_roots]
    prev_roots_diag_bound = [-k*vd.diagonal_values[1]*prev_root[1] for prev_root in prev_roots_diag]
    prev_roots_diag_last_non_zero_coordinate = [findlast(≠(0),prev_root) for prev_root in prev_roots_diag]
    prev_roots_constraints = (prev_roots_diag,prev_roots_diag_bound,prev_roots_diag_last_non_zero_coordinate)

    #@info "roots_for_pair($pair,$prev_roots)"
    roots = extend_root_stem(vd,[k],k .* vd.diagonal_basis[1],l,l-k^2*vd.diagonal_values[1],prev_roots_constraints,t2_cache)
    
    @assert all(times(vd,root,prev) ≤ 0 for root in roots for prev in prev_roots)  "All angles with previous roots should be acute."
    # just in case
    filter!(root -> all(times(vd,root,prev) ≤ 0 for prev in prev_roots),roots)
    @assert all(is_root(vd.quad_space,vd.ring,root) for root in roots) "All outputs of extend_root_stem must be roots"
    @assert all(norm_squared(vd,root) == l for root in roots) "All outputs of extend_root_stem must have correct length"
    @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
    @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂) """
    Two roots at same distance to basepoint, and both acute with previous ones, should be acute with each other:
    * why?
    * distance is $fake_dist ≈ $(approx(fake_dist)))
    * roots are  $roots
    * gram = $([approx(Gram_coeff(vd.quad_space,r₁,r₂)) for r₁ in roots for r₂ in roots if r₁≠r₂])
    * angles = $([(times(vd,r₁,r₂),approx(times(vd,r₁,r₂))) for r₁ in roots for r₂ in roots if r₁≠r₂])
    """
   
    
    return roots
    
end

function roots_for_next_pair!(vd,dict,prev_roots;t2_cache=nothing)

    pair = next_min_pair!(vd,dict)
    k,l = pair
    fake_dist = -(k^2)*vd.diagonal_values[1].elem_in_nf//(l.elem_in_nf)
    @info "next pair is $pair (fake_dist is $fake_dist ≈ $(approx(fake_dist)))"


    # Here we should assert that the roots in prev_roots are actually closer than the next min pair can give us, otherwise we will get inconsistencies

    roots =  roots_for_pair(vd,pair,prev_roots,t2_cache=t2_cache)

    @toggled_assert all(fake_dist_to_basepoint(vd,r) == fake_dist for r in roots)
    return roots

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

        Cox_matrix = Coxeter_matrix(vd.quad_space, vd.ring, prev_roots) 
        das = build_diagram_and_subs(Cox_matrix,vd.dim-1)
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
