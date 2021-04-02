
mutable struct VinbergData

    dim::Int
    field::AnticNumberField
    ring::NfAbsOrd{AnticNumberField,nf_elem}
    gram_matrix::AbstractAlgebra.Generic.MatSpaceElem{nf_elem}
    quad_space::Hecke.QuadSpace{AnticNumberField,AbstractAlgebra.Generic.MatSpaceElem{nf_elem}}

    diagonal_basis::Vector{Vector{nf_elem}}
    diagonal_values::Vector{nf_elem}
    scaling::Vector{nf_elem}
    
    diagonal_change::Array{nf_elem,2}
    diagonal_change_inv::Array{nf_elem,2}


    possible_root_norms_squared_up_to_squared_units::Vector{nf_elem}

    # Precomputed stuff:
    diago_over_scaling::Vector{nf_elem}
    diago_over_scalingsq::Vector{nf_elem}
    two_diago_over_scaling_times_length::Dict{nf_elem,Vector{nf_elem}}
    diago_vector_last_on_coordinates::Vector{Vector{Int}}

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

    lengths = possible_root_norms_squared_up_to_squared_units(ring_of_integers, number_field, quad_space) .|> (x -> x.elem_in_nf)
   
    # Precomputations ---
    diago_over_scaling = [α.elem_in_nf//s.elem_in_nf for (α,s) in zip(diagonal_values,scaling)]::Vector{nf_elem} 
    diago_over_scalingsq = [α.elem_in_nf//(s.elem_in_nf^2) for (α,s) in zip(diagonal_values,scaling)]::Vector{nf_elem} 
    two_diago_over_scaling_times_length = Dict([l => (2//l) .* diago_over_scaling for l in lengths])::Dict{nf_elem,Vector{nf_elem}}
    
    is_last_on_coord(vecs,v_i,c_i) = vecs[v_i][c_i]≠0 && all(vecs[i][c_i] == 0 for i in v_i+1:length(vecs))
    diago_vector_last_on_coordinates =[filter(c_i -> is_last_on_coord(diagonal_basis_vecs,v_i,c_i), 1:n) for v_i in 1:n] ::Vector{Vector{Int}}
    # -------------------

    vd =  VinbergData(
        n,
        number_field,
        ring_of_integers,
        matrix(number_field,gram_matrix),
        quad_space,
        [[v.elem_in_nf for v in vec] for vec in diagonal_basis_vecs],
        [v.elem_in_nf for v in diagonal_values],
        [v.elem_in_nf for v in scaling],
        Matrix(matrix(number_field,diagonal_basis_vecs))',    # notice the transpose here and below. Much frustration took place before I found out I needed those!
        Matrix(inv(matrix(number_field,diagonal_basis_vecs)))',
        lengths,
        diago_over_scaling,
        diago_over_scalingsq,
        two_diago_over_scaling_times_length,
        diago_vector_last_on_coordinates,
    )

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

is_root(vd::VinbergData, vector::Vector) = is_root(vd.quad_space,vd.ring,vector)
Gram_coeff(vd::VinbergData,r₁::Vector,r₂::Vector) = Gram_coeff(vd.quad_space,r₁,r₂)
Coxeter_coeff(vd::VinbergData,r₁::Vector,r₂::Vector) = Coxeter_coeff(vd.quad_space,vd.ring,r₁,r₂)
Coxeter_matrix(vd,roots) = Coxeter_matrix(vd.quad_space,vd.ring,roots)

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
    nf_elem,
    Tuple{
        nf_elem,
        Vector{nf_elem}
    }
}



function enumerate_k(vd::VinbergData,l,k_min,k_max)
    
    @info "enumerate_k (l=$l, k_min=$k_min, k_max=$k_max)"

    ring = vd.ring
    field = vd.field

    α_1 = vd.diagonal_values[1]
    s_1 = vd.scaling[1]

    @assert α_1 < 0 

    M = approx_sum_at_places(field(l*s_1^2)//field(α_1),first_place_idx=2)
    P = infinite_places(field)
    k_min_squared_approx = approx(field((k_min*s_1)^2))
    k_max_squared_approx = approx(field((k_max*s_1)^2))

    upscaled_candidates = non_neg_short_t2_elems(ring, k_min_squared_approx-1 , k_max_squared_approx  + M + 1)
    
    if upscaled_candidates[end][1] == 0
        pop!(upscaled_candidates)
    end

    candidates::Vector{nf_elem} = map(x -> x[1]//s_1, upscaled_candidates) # Correct the scaling  
    # We only need k>0 because k=0 corresponds to hyperplane containing the basepoint, treated elsewhere
    @toggled_assert all(>(0),candidates)

    
   
    # crystallographic_condition
    filter!(
        k -> divides(l,2*k*α_1,ring),
        candidates,    
    )
    
    # conjugates are of bounded length
    filter!(
        k -> all(≤(α_1*k^2,l,p) for p in P[2:end]),
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
            l => next_min_k_for_l(vd,0,[],l)
            for l in vd.possible_root_norms_squared_up_to_squared_units
        ]
    )
end


function next_min_pair(vd,dict)
    field = vd.field
    min_pair = nothing

    # that's the value we want to minimize (related to the distance)
    val(k,l) = k^2//l
    

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
    elems::Vector{Vector{Tuple{nf_elem,fmpq,Vector{Float64}}}}
end

function BoundedT2ElemsCache(field)
    return (BoundedT2ElemsCache)(Vector{Int}([0]),Vector{Vector{Tuple{nf_elem,fmpq,Vector{Float64}}}}([[(field(0),fmpq(0),Float64.(conjugates_real(field(0))))]]))
end


function bounded_t2_elems!(
    field,
    ring,
    t2_bound,
    cache,
)
    int_bound = ceil(t2_bound)
    if int_bound > cache.bounds[end]
        with_margin = ceil(int_bound*1.2)
        push!(cache.bounds,with_margin)
        new_elems = filter(
            x -> x∉cache.elems[end],
            non_neg_short_t2_elems(ring,cache.bounds[end-1],cache.bounds[end]),# .|> abs,
        )
        sort!(new_elems,by=(x->x[1]))
        new_elems_nf = map(x -> (x[1].elem_in_nf,x[2],Float64.(conjugates_real(x[1].elem_in_nf))), new_elems)::Vector{Tuple{nf_elem,fmpq,Vector{Float64}}}
        push!(cache.elems, new_elems_nf)
    end
    
    return searchsortedfirst(cache.bounds,t2_bound)
end



const AffineConstraints = Tuple{
    Vector{Vector{nf_elem}}, # Previous roots
    Vector{nf_elem},         # ≤ bound
    Vector{Int},           # last non non-zero coordinate
}     


function MkAffineConstraints(
    vecs::Vector{Vector{nf_elem}}, # Previous roots
    vals::Vector{nf_elem},         # ≤ bound
)
    @assert all(any(≠(0),vec) for vec in vecs)
    last_nz = Int[findlast(≠(0),vec) for vec in vecs]
    return (vecs,vals,last_nz)
end

no_constraints() = MkAffineConstraints(Vector{Vector{nf_elem}}(),Vector{nf_elem}())

function clearly_inconsistent(ac::AffineConstraints,idx,dim)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac 
    return any( bound < 0 && last_non_zero < idx for (bound,last_non_zero) in zip(c_values,c_last_non_zero_coordinates))   
end

function update_constraints(ac::AffineConstraints,idx::Int,val_at_idx::nf_elem)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac
    
    return (c_vectors,[b - val_at_idx*r[idx] for (r,b) in zip(c_vectors,c_values)],c_last_non_zero_coordinates)

end


const Interval = Tuple{Tuple{Bool,nf_elem},Tuple{Bool,nf_elem}}

function interval_for_k_j(field,ac::AffineConstraints,j::Int)::Interval
    applicable = [(vec[j],val) for (vec,val,last) in zip(ac...) if last == j]
    pos = [v//c for (c,v) in applicable if c > 0]
    neg = [v//c for (c,v) in applicable if c < 0]
    no_lb = isempty(neg)
    no_ub = isempty(pos)
    lb = ( no_lb ? field(1) : maximum(neg) )::nf_elem
    ub = ( no_ub ? field(0) : minimum(pos) )::nf_elem
    return ((no_lb,lb),(no_ub,ub))
end

# does not reduce allocations.
# We would need to use AbstractAlgebra's dangerous mul! and others I think
function interval_for_k_j2(field,ac::AffineConstraints,j::Int)::Interval
    vecs,vals,last = ac
    @assert length(vecs) == length(vals)
    @assert length(vecs) == length(last)

    lb = field(1)
    ub = field(0)
    no_lb = true
    no_ub = true
    @inbounds for idx in 1:length(vecs)
        if last[idx] == j 
            c = vecs[idx][j]
            v = vals[idx]
            a = v//c
            # div!(a,v,c) this function does not exist
            if c > 0
                no_ub && (ub = a)
                !no_ub && (ub = min(ub,a))
                no_ub = false
            elseif c < 0
                no_lb && (lb = a)
                !no_lb && (lb = max(lb,a))
                no_lb = false
            else
                #@assert false "unreachable"
            end
        end
    end
    return ((no_lb,lb),(no_ub,ub))
end

function in_interval(
    k,
    interval::Interval
)
    ((no_lb,lb),(no_ub,ub)) = interval
    return (no_lb || k ≥ lb) && (no_ub || k ≤ ub)
end

function is_empty(i::Interval)
    (no_lb,lb),(no_ub,ub) = i
    return !no_lb && !no_ub && (lb > ub)
end

@inline function _extend_root_stem_full!(
    vd::VinbergData,
    stem::Vector{nf_elem},
    stem_can_rep::Vector{nf_elem},
    stem_length::Int,
    root_length::nf_elem,
    root_length_minus_stem_norm_squared::nf_elem,
    constraints::AffineConstraints,
    t2_cache::BoundedT2ElemsCache,
    roots::Vector{Vector{nf_elem}}
)

    field = vd.field
    ring = vd.ring
    space = vd.quad_space
    P = infinite_places(field)
    l = root_length   
    (c_vectors,c_values,c_last_non_zero_coordinates) = constraints

    @toggled_assert stem_can_rep == to_can_rep(vd,stem) "Sanity check. `stem_can_rep` must always be equal to `to_can_rep(vd,stem)`!"
    @toggled_assert times(vd,stem_can_rep,stem_can_rep) == l "Sanity check. If we are here, by the previous case (j==vd.dim) necessarily we have the right length."
    @toggled_assert all(bound ≥ 0 for bound in c_values) "All affine constraints given by previous roots must be satisfied."

    if is_integral(space, ring, stem_can_rep) && is_root(space,ring,stem_can_rep,l) 
        # integralness should be guaranteed to hold for all but the last coordinate I think.
        append!(roots,[deepcopy(stem_can_rep)])
    end
end

@inline function _extend_root_stem_one_coord_left!(
    vd::VinbergData,
    stem::Vector{nf_elem},
    stem_can_rep::Vector{nf_elem},
    stem_length::Int,
    root_length::nf_elem,
    root_length_minus_stem_norm_squared::nf_elem,
    constraints::AffineConstraints,
    t2_cache::BoundedT2ElemsCache,
    roots::Vector{Vector{nf_elem}}
)
   
    
    field = vd.field
    ring = vd.ring
    space = vd.quad_space
    P = infinite_places(field)
    l = root_length

    j = stem_length + 1

    α_j = vd.diagonal_values[j]
    v_j = vd.diagonal_basis[j]    
    s_j = vd.scaling[j]
    l_j = root_length_minus_stem_norm_squared # l - sum([vd.diagonal_values[i]*stem[i]^2 for i in 1:length(stem)]) 
    
    #α_over_s = α_j//s_j
    α_over_s = vd.diago_over_scaling[j]
    #α_over_s² = α//s_j^2
    α_over_s² = vd.diago_over_scalingsq[j]
    #two_α_over_ls = 2*α//(l*s_j)
    two_α_over_ls = vd.two_diago_over_scaling_times_length[l][j]



    # If k₁,…,k_{j-1} are chosen and k_j is the last one that needs to be found.
    # If (α₁,…,α_j) is the diagonalized inner product, and r = (k₁,…,k_j) the root to be found, of normed² == l, then
    # We have
    # 
    #   ∑_{i=1}^{j-1} k_i^2α_i + k_j^2α_j == normed²(r) == l 
    #
    # Which means k_j^2 = (l - ∑_{i=1}^{j-1}k_i^2α_i)/α_j.
    issquare,square_root = Hecke.issquare(l_j//α_j)
    
    
    if issquare && divides(l,2*square_root*α_j,ring) # crystal
        
        k = square_root
        stem_updated = deepcopy(stem) # This is very important for correctness
        stem_can_rep_updated = deepcopy(stem_can_rep) # Same

        stem_updated[j] = k
        #stem_can_rep_updated = stem_can_rep .+ k .* v_j
        u_plus_k_v(stem_can_rep_updated,stem_can_rep,k,v_j)

        _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,l,l_j - k^2*α_j,update_constraints(constraints,j,k*α_j),t2_cache,roots)
        if k ≠ 0
            stem_updated[j] = -k
            #stem_can_rep_updated = stem_can_rep .- k .* v_j
            u_plus_k_v(stem_can_rep_updated,stem_can_rep,k,-v_j)
            _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,l,l_j - k^2*α_j,update_constraints(constraints,j,-k*α_j),t2_cache,roots)
        end
        
    else
        return 
    end
   
end


@inline function find_range(
    field,
    interval_αk::Interval,
    α_over_s::nf_elem,
    ordered::Vector{Tuple{nf_elem,fmpq,Vector{Float64}}}
   )::Tuple{Int,Int,Int,Int}
   
    dummyvec = Float64[]

    (no_lb,lb),(no_ub,ub) = interval_αk
    lb  = interval_αk[1][2]//α_over_s 
    ub  = interval_αk[2][2]//α_over_s 
    
    if (no_lb || lb ≤ 0) && (no_ub || ub ≥ 0) 

        last_idx_neg = (no_lb ? length(ordered) : searchsortedlast(ordered,(-lb,:dummy,dummyvec),by=(x->x[1])))
        last_idx_pos = (no_ub ? length(ordered) : searchsortedlast(ordered,(ub,:dummy,dummyvec),by=(x->x[1])))
        return (1,last_idx_pos,1,last_idx_neg)

    elseif (!no_lb && lb ≥ 0) && (no_ub || ub ≥ 0)
        
        first_idx_lb = (no_lb ? 1 : searchsortedfirst(ordered,(lb,:dummy,dummyvec),by=(x->x[1])))
        last_idx_ub = (no_ub ? length(ordered) : searchsortedlast(ordered,(ub,:dummy,dummyvec),by=(x->x[1])))
        return (first_idx_lb,last_idx_ub,1,0)


    elseif (no_lb || lb ≤ 0) && (!no_ub && ub ≤ 0)
        
        first_idx_ub = (no_ub ? 1 : searchsortedfirst(ordered,(-ub,:dummy,dummyvec),by=(x->x[1])))
        last_idx_lb = (no_lb ? length(ordered) : searchsortedlast(ordered,(-lb,:dummy,dummyvec),by=(x->x[1])))
        return (1,0,first_idx_ub,last_idx_lb)
    end
    
    @assert false "should not be reachable"

end

function _extend_root_stem!(
    vd::VinbergData,
    stem::Vector{nf_elem},
    stem_can_rep::Vector{nf_elem},
    stem_length::Int,
    root_length::nf_elem,
    root_length_minus_stem_norm_squared::nf_elem,
    constraints::AffineConstraints,
    t2_cache::BoundedT2ElemsCache,
    roots::Vector{Vector{nf_elem}}
)
 
    
    j = stem_length + 1 
    
    if clearly_inconsistent(constraints,j,vd.dim)
        return
    end
    

    if j == vd.dim + 1
        return _extend_root_stem_full!(vd,stem,stem_can_rep,stem_length,root_length,root_length_minus_stem_norm_squared,constraints,t2_cache,roots)
    end

    if j == vd.dim
        return _extend_root_stem_one_coord_left!(vd,stem,stem_can_rep,stem_length,root_length,root_length_minus_stem_norm_squared,constraints,t2_cache,roots)
    end

    field = vd.field
    ring = vd.ring
    space = vd.quad_space
    P = infinite_places(field)
    l = root_length
    (c_vectors,c_values,c_last_non_zero_coordinates) = constraints

    α_j = vd.diagonal_values[j]
    v_j = vd.diagonal_basis[j]    
    s_j = vd.scaling[j]
    l_j = root_length_minus_stem_norm_squared # l - sum([vd.diagonal_values[i]*stem[i]^2 for i in 1:length(stem)]) 
    
    #α_over_s = α_j//s_j
    α_over_s = vd.diago_over_scaling[j]
    #α_over_s² = α//s_j^2
    α_over_s² = vd.diago_over_scalingsq[j]
    #two_α_over_ls = 2*α//(l*s_j)
    two_α_over_ls = vd.two_diago_over_scaling_times_length[l][j]


  
    
    #t2_bound_for_sk = approx_sum_at_places(l_j//(α_over_s²),first_place_idx=1)+1
    t2_bound_for_sk = Hecke.tr(l_j//(α_over_s²))
    last_bounded_t2_candidates_vector_idx = bounded_t2_elems!(
        vd.field,
        vd.ring, 
        t2_bound_for_sk,
        t2_cache
    )
   
    conjugates_bound = Float64.(conjugates_real(l_j//α_over_s²))
   
    #    Constraint on norm at al places:
    #
    #       k^2α ≤ l_j at all places
    #    ⇔  (sk)^2 α/s^2 ≤ l_j at all places
    #    ⇔  (sk)^2 * α_over_s² ≤ l_j at all places
    #
    #good_norm(sk,conjs_sk) = all( ≤(sk^2,l_j // α_over_s²,p) for p in P)
    good_norm(sk,conjs_sk) = all( conjs_sk[i]^2 ≤ conjugates_bound[i]+1 for i in 1:length(conjugates_bound)) # that's approximative but should be good enough
    #good_norm(sk) = true # all( ≤(sk^2 * α_over_s²,l_j,p) for p in P)

    #    Constraint as given by the crystallographic condition:
    #
    #       l | 2kα  
    #    ⇔  2kα/l ∈ ring  
    #    ⇔  sk (2α//ls) ∈ ring 
    #    ⇔  sk (two_α_over_ls) ∈ ring 
    crystal(sk) = sk * two_α_over_ls ∈ ring

    integral(a_stem_can_rep) = all((a_stem_can_rep[idx] ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j])

    
    # The idea is that interval_k_j gives an interval outside of which k_jα_j is not valid due to the constraints of acute angles given by previous roots.
    # The code below SHOULD then use this interval to only iterate over k_js in this interval.
    interval_αk = interval_for_k_j(field,constraints,j) 
    

    # If the endpoints are not ±∞, rescale them to get endpoints for sk instead of endpoints of α*k
    @assert α_j > 0
    if is_empty(interval_αk)
        return 
    end

    stem_updated = deepcopy(stem)
    stem_can_rep_updated = deepcopy(stem_can_rep) #copy(stem_can_rep)

    
    for (idx,ordered) in enumerate(t2_cache.elems[1:last_bounded_t2_candidates_vector_idx])
       
        @toggled_assert issorted(ordered)
        isempty(ordered) && continue
        
        (first_idx_pos,last_idx_pos,first_idx_neg,last_idx_neg) = find_range(field,interval_αk, α_over_s, ordered) 
        for i in min(first_idx_pos,first_idx_neg):max(last_idx_pos,last_idx_neg)
            
            sk,t2sk,conjs_sk = ordered[i]
            pos =  first_idx_pos ≤ i && last_idx_pos ≥ i
            neg = first_idx_neg ≤ i && last_idx_neg ≥ i
            check_t2 = idx == last_bounded_t2_candidates_vector_idx

            if (check_t2 ⇒ (t2sk ≤ t2_bound_for_sk)) &&
                crystal(sk) &&
                good_norm(sk,conjs_sk)
                #sk * two_α_over_ls ∈ ring && # crystallographic condition 
                #all( ≤(sk^2 * α_over_s²,l_j,p) for p in P) #  norms are OK
                
                k = sk // s_j
                if pos 
                    @toggled_assert in_interval(k*α_j,interval_αk)
                    stem_updated = copy!(stem_updated,stem); stem_updated[j] = k
                    #stem_can_rep_updated = stem_can_rep .+ k .* v_j
                    u_plus_k_v(stem_can_rep_updated,stem_can_rep,k,v_j)
                    if integral(stem_can_rep_updated) # all((stem_can_rep_updated[idx] ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j]) # integral
                        _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,l,l_j - k^2*α_j,update_constraints(constraints,j,k*α_j),t2_cache,roots)
                    end
                end
                if neg && k≠0 
                    @toggled_assert in_interval(-k*α_j,interval_αk)
                    stem_updated = copy!(stem_updated,stem); stem_updated[j] = -k
                    #stem_can_rep_updated = stem_can_rep .- k .* v_j
                    u_plus_k_v(stem_can_rep_updated,stem_can_rep,-k,v_j)
                    if integral(stem_can_rep_updated) # all((stem_can_rep_updated[idx] ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j]) # integral
                        _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,l,l_j - k^2*α_j,update_constraints(constraints,j,-k*α_j),t2_cache,roots)
                    end
                end

            end

        end
    end



end


function u_plus_k_v(to,u,k,v)
    @assert length(to) == length(u)
    @assert length(to) == length(v)
    @inbounds for i in 1:length(to)
        mul!(to[i],k,v[i])
        addeq!(to[i],u[i])
    end
    @toggled_assert to == u + (k .* v)
end

function extend_root_stem(
    vd::VinbergData,
    stem_diag_rep::Vector{nf_elem},
    root_length,
    constraints::AffineConstraints,
    t2_cache::BoundedT2ElemsCache
)

    
    stem_length = length(stem_diag_rep)
    stem_diag_rep = vcat(stem_diag_rep,fill(vd.field(0),vd.dim-length(stem_diag_rep)))
    stem_can_rep = to_can_rep(vd,stem_diag_rep)
    
    stem_norm_squared = norm_squared(vd,stem_can_rep)
    
    j = 1
    #lol#println("| "^(j-1))
    #lol#println("| "^(j-1), "---------------------------") 
    #lol#println("| "^(j-1), "extend_root_stem")
    #lol#println("| "^(j-1), "stem_diag_rep             :  ", stem_diag_rep)
    #lol#println("| "^(j-1), "of length                 :  ", stem_length)
    #lol#println("| "^(j-1), "stem_can_rep              :  ", stem_can_rep)
    #lol#println("| "^(j-1), "root length               :  ", root_length)
    #lol#println("| "^(j-1), "root length_minu_partial  :  ", root_length - stem_norm_squared)
    #lol#println("| "^(j-1), "constraints               :  ", constraints)
    #lol#println("| "^(j-1), "---------------------------") 



    #@info "roots_for_pair($pair,$prev_roots)"
    roots_go_here = Vector{Vector{nf_elem}}()
     _extend_root_stem!(vd,stem_diag_rep,stem_can_rep,stem_length,root_length,root_length-stem_norm_squared,constraints,t2_cache,roots_go_here)
     return roots_go_here
end

function roots_at_distance_zero(vd::VinbergData)

    t2_cache = BoundedT2ElemsCache(vd.field) 
    zero_stem = nf_elem[vd.field(0)]
    
    return vcat([extend_root_stem(vd,zero_stem,l,no_constraints(),t2_cache) for l in vd.possible_root_norms_squared_up_to_squared_units]...)
end

function cone_roots(vd,roots_at_distance_zero)

    @warn "Cone roots computation are approximative ⇒ double check the results by hand."
    @warn "If results look dubious, increasing LP precision by setting VinbergsAlgorithmNF.LP_PRECISION to something higher (currently $LP_PRECISION)"
    roots_at_distance_zero = [root for root in roots_at_distance_zero]
    @info "starting with $(length(roots_at_distance_zero)) roots at distance zero"
    
    len = length(roots_at_distance_zero)

    # We put first the roots with integer coordinates to maximize the chance of having them in the cone roots
    # It's not necessary but easier to analyze the output and compare with rgug then
    integer_roots = filter(r -> all(isinteger,r), roots_at_distance_zero)
    non_integer_roots = filter(r -> !all(isinteger,r), roots_at_distance_zero)
    sort!(integer_roots)
    roots_at_distance_zero = vcat(integer_roots,non_integer_roots)
    @assert len == length(roots_at_distance_zero) "did we drop a root while reorgarizing them?"
    
    cone_roots = Vector{Vector{nf_elem}}()


    for r in roots_at_distance_zero
        @debug "looking at $r"
        if  all((-1)*r ≠ cr for cr in cone_roots)
            @debug "so far so good"
            if is_necessary_halfspace(vd.gram_matrix.entries,cone_roots,-r)
                @debug "degeneration"
                push!(cone_roots,r)
            end
        
        end
        @debug "have $(length(cone_roots)) cone roots" 
    end
    
    cone_roots = drop_redundant_halfspaces(vd.gram_matrix.entries,cone_roots)
    @assert all(times(vd,r₁,r₂) ≤ 0 for r₁ in cone_roots for r₂ in cone_roots if r₁≠r₂)

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
        t2_cache = BoundedT2ElemsCache(vd.field) 
    end

    (k,l) = pair
    fake_dist = -(k^2)*vd.diagonal_values[1]//(l)
   
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
    prev_roots_constraints = MkAffineConstraints(prev_roots_diag,prev_roots_diag_bound)

    #@info "roots_for_pair($pair,$prev_roots)"
    roots = extend_root_stem(vd,[k],l,prev_roots_constraints,t2_cache)
    
    @assert all(times(vd,root,prev) ≤ 0 for root in roots for prev in prev_roots)  "All angles with previous roots should be acute."
    @assert all(is_root(vd,root) for root in roots) "All outputs of extend_root_stem must be roots"
    @assert all(norm_squared(vd,root) == l for root in roots) "All outputs of extend_root_stem must have correct length"
    @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
    @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂) "Two roots at same distance and compatible with everything before them should be compatible with each other." 
   
    
    return roots
    
end

function roots_for_next_pair!(vd,dict,prev_roots;t2_cache=nothing)

    pair = next_min_pair!(vd,dict)
    k,l = pair
    fake_dist = -(k^2)*vd.diagonal_values[1]//l
    @info "next pair is $pair (fake_dist is $fake_dist ≈ $(approx(fake_dist,8)))"


    # Here we should assert that the roots in prev_roots are actually closer than the next min pair can give us, otherwise we will get inconsistencies

    roots =  roots_for_pair(vd,pair,prev_roots,t2_cache=t2_cache)

    @toggled_assert all(fake_dist_to_basepoint(vd,r) == fake_dist for r in roots)
    return roots

end

function next_n_roots!(vd,prev_roots,dict,das;n=10,t2_cache=nothing)

    if t2_cache===nothing
        t2_cache = BoundedT2ElemsCache(vd.field)
    end
    
    # TODO
    # asserts that the previous roots are closer than what the dict proposes, otherwise we'll get inconsistencies


    roots = prev_roots
    #Coxeter_matrix = get_Coxeter_matrix(vd.quad_space, vd.ring, prev_roots) 
    new_roots = Vector{Vector{nf_elem}}()
    while n > 0 


        new_roots = roots_for_next_pair!(vd,dict,roots;t2_cache=t2_cache)
        n = n - length(new_roots)

        for root in new_roots 
            @info "Got new root $root"
            extend!(das,[Coxeter_coeff(vd.quad_space, vd.ring, old,root) for old in roots])
            push!(roots,root)
        end


        if is_finite_volume(das)
            # Sanity checks
            @assert all(is_root(vd.quad_space,vd.ring,r) for r in roots)
            @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂)
            @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
            @assert all(fake_dist_to_basepoint(vd,roots[i]) ≤ fake_dist_to_basepoint(vd,roots[i+1]) for i in 1:length(roots)-1)
            return (true,(roots,dict,das))
        end
    end

    # Sanity checks
    @assert all(is_root(vd.quad_space,vd.ring,r) for r in roots)
    @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂)
    @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
    @assert all(fake_dist_to_basepoint(vd,roots[i]) ≤ fake_dist_to_basepoint(vd,roots[i+1]) for i in 1:length(roots)-1)

    return (false,(roots,dict,das))
end

function next_n_roots!(
    vd::VinbergData,
    prev_roots::Vector;
    n=10,
    t2_cache=nothing
)

    # TODO: make this work even when the previous roots are not cone roots.
    # In this case, one just need to iterate over the pairs (k,l) in the dictionary until we're at distance ≥ than the farthest root in prev_roots
    @assert all(fake_dist_to_basepoint(vd,r) == 0 for r in prev_roots) "Only works with previous roots == cone roots"

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
