
struct VinbergData

    dim::Int
    field::AnticNumberField
    ring::NfAbsOrd{AnticNumberField,nf_elem}
    gram_matrix::AbstractAlgebra.Generic.MatSpaceElem{nf_elem}
    quad_space::Hecke.QuadSpace{AnticNumberField,AbstractAlgebra.Generic.MatSpaceElem{nf_elem}}

    diagonal_basis::Vector{Vector{boxed_nf_elem}}
    diagonal_values::Vector{boxed_nf_elem}
    scaling::Vector{boxed_nf_elem}
    
    diagonal_change::Array{nf_elem,2}
    diagonal_change_inv::Array{nf_elem,2}


    possible_root_norms_squared_up_to_squared_units::Vector{boxed_nf_elem}

    # Precomputed stuff:
    diago_over_scaling::Vector{boxed_nf_elem}
    diago_over_scalingsq::Vector{boxed_nf_elem}
    two_diago_over_scaling_times_length::Dict{Int,Vector{boxed_nf_elem}}
    diago_vector_last_on_coordinates::Vector{Vector{Int}}

end


function VinbergData(number_field,gram_matrix,v₀=nothing)

    (n,m) = size(gram_matrix)
    @assert n == m "The Gram gram_matrix must be square."

    ring_of_integers = maximal_order(number_field)
    quad_space = quadratic_space(number_field, matrix(number_field,gram_matrix))
    #quad_lattice = Hecke.lattice(quad_space)
        
    @assert is_feasible(quad_space) "The quadratic form must be feasible sig (n,1) and all conjugates sig (n+1,0)"
    @assert all(number_field(c) ∈ ring_of_integers for c in gram_matrix) "The Gram matrix must have coefficients in the ring of integers."

    if v₀≠nothing
        @assert v₀' * gram_matrix * v₀ < 0
    end

    diagonal_basis_vecs,diagonal_values,scaling = diagonalize_and_get_scaling(gram_matrix,ring_of_integers,number_field,v₀)
    negative_vector_index = filter(x-> x[2]<0, collect(enumerate(diagonal_values)))[1][1]
    
    if negative_vector_index ≠ 1
        diagonal_basis_vecs[1],diagonal_basis_vecs[negative_vector_index] = diagonal_basis_vecs[negative_vector_index],diagonal_basis_vecs[1]
        diagonal_values[1],diagonal_values[negative_vector_index] = diagonal_values[negative_vector_index],diagonal_values[1]
        scaling[1],scaling[negative_vector_index] = scaling[negative_vector_index],scaling[1]
    end

    if v₀≠nothing
        @assert v₀ == diagonal_basis_vecs[1]
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
    two_diago_over_scaling_times_length = Dict([idx => box.((2//l) .* diago_over_scaling) for (idx,l) in enumerate(lengths)])::Dict{Int,Vector{boxed_nf_elem}}
    
    is_last_on_coord(vecs,v_i,c_i) = vecs[v_i][c_i]≠0 && all(vecs[i][c_i] == 0 for i in v_i+1:length(vecs))
    diago_vector_last_on_coordinates =[filter(c_i -> is_last_on_coord(diagonal_basis_vecs,v_i,c_i), 1:n) for v_i in 1:n] ::Vector{Vector{Int}}
    # -------------------

    vd =  VinbergData(
        n,
        number_field,
        ring_of_integers,
        matrix(number_field,gram_matrix),
        quad_space,
        [[box(v.elem_in_nf) for v in vec] for vec in diagonal_basis_vecs],
        [box(v.elem_in_nf) for v in diagonal_values],
        [box(v.elem_in_nf) for v in scaling],
        Matrix(matrix(number_field,diagonal_basis_vecs))',    # notice the transpose here and below. Much frustration took place before I found out I needed those!
        Matrix(inv(matrix(number_field,diagonal_basis_vecs)))',
        box.(lengths),
        box.(diago_over_scaling),
        box.(diago_over_scalingsq),
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
    @assert ex(vd.diagonal_values[1]) < 0 "sanity check"

    return ex.(vd.diagonal_basis[1])
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
    @toggled_assert fake_dist == -(to_diag_rep(vd,u)[1]^2)*ex(vd.diagonal_values[1])//norm_squared(vd,u) 

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
    
    @debug "enumerate_k (l=$l, k_min=$k_min, k_max=$k_max)"

    ring = vd.ring
    field = vd.field

    α_1 = ex(vd.diagonal_values[1])
    s_1 = ex(vd.scaling[1])

    @assert α_1 < 0 

    M = approx_sum_at_places(field(ex(l)*s_1^2)//field(α_1),first_place_idx=2)
    P = infinite_places(field)
    k_min_squared_approx = approx(field((k_min*s_1)^2))
    k_max_squared_approx = approx(field((k_max*s_1)^2))

    upscaled_candidates = non_neg_short_t2_elems(ring, k_min_squared_approx-1 , k_max_squared_approx  + M + 1)
    
    if upscaled_candidates[end][1] == 0
        pop!(upscaled_candidates)
    end

    candidates::Vector{nf_elem} = map(x -> x[1]//s_1, upscaled_candidates) # Correct the scaling  
    # We only need k>0 because k=0 corresponds to hyperplane containing the basepoint, treated elsewhere
    @toggled_assert all(x->x>0,candidates)

    
   
    # crystallographic_condition
    filter!(
            k -> divides(ex(l),2*k*α_1,ring),
        candidates,    
    )
    
    # conjugates are of bounded length
    filter!(
            k -> all(≤(α_1*k^2,ex(l),p) for p in P[2:end]),
        candidates,
    )
      
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
    val(k,l) = k^2//ex(l)
    

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
    lock::ReentrantLock
    bounds::Vector{Int}
    elems::Vector{Vector{Tuple{
                               boxed_nf_elem,
                               Vector{Interval{Float64}},
                               fmpq,
                               BitMatrix
                              }
                        }
                 }
end

function BoundedT2ElemsCache(vd::VinbergData)
    field = vd.field
    return BoundedT2ElemsCache(ReentrantLock(),Vector{Int}([0]),Vector{Vector{Tuple{boxed_nf_elem,Vector{Interval{Float64}},fmpq,BitMatrix}}}([[(box(field(0)),get_enclosing_intervals_float64(field(0)),fmpq(0),crystal_matrix(vd,field(0)))]]))
end

function crystal_matrix(vd::VinbergData,sk)
    return BitMatrix(sk* ex(vd.two_diago_over_scaling_times_length[l][j]) ∈ vd.ring for l in 1:length(vd.possible_root_norms_squared_up_to_squared_units), j in 1:vd.dim)
end

function bounded_t2_elems!(
    vd::VinbergData,
    t2_bound,
    cache,
)
    ring = vd.ring
    field = vd.field

    int_bound = ceil(t2_bound)
    if int_bound > cache.bounds[end]
        lock(cache.lock) do
            with_margin = ceil(int_bound*1.2)
            push!(cache.bounds,with_margin)
            new_elems = filter( # could probably me optimized to not have to filter away
                x -> all(x[2] ≠ cache.bounds[end-1]),
                non_neg_short_t2_elems(ring,cache.bounds[end-1],cache.bounds[end]),# .|> abs,
            )
            @toggled_assert all(x≥0 for (x,t) in new_elems)
            @toggled_assert all(all(x[1].elem_in_nf ≠ y[1] for y in cache.elems[end]) for x in new_elems)
            
            # for each elem sk, we keep sk, its T₂ norm, lower bounds on the absolute values of the conjugates of sk, and each possible divisibility condition to check the crystallographic condition
            new_elems_nf = map(x -> (box(x[1].elem_in_nf),get_enclosing_intervals_float64(x[1].elem_in_nf),x[2],crystal_matrix(vd,x[1].elem_in_nf)), new_elems)
            sort!(new_elems_nf,by=(x->lo(x[1])))
            sort!(new_elems_nf,by=(x->ex(x[1])))
            push!(cache.elems, new_elems_nf)
        end
    end
    
    return searchsortedfirst(cache.bounds,t2_bound)
end



const AffineConstraints = Tuple{
    Vector{Vector{boxed_nf_elem}}, # Previous roots
    Vector{boxed_nf_elem},         # ≤ bound
    Vector{Int},           # last non non-zero coordinate
}     


function MkAffineConstraints(
    vecs::Vector{Vector{boxed_nf_elem}}, # Previous roots
    vals::Vector{boxed_nf_elem},         # ≤ bound
)
    @assert all(any(x->ex(x)≠0,vec) for vec in vecs)
    last_nz = Int[findlast(x->ex(x)≠0,vec) for vec in vecs]
    return (vecs,vals,last_nz)
end

no_constraints() = MkAffineConstraints(Vector{Vector{boxed_nf_elem}}(),Vector{boxed_nf_elem}())

function clearly_inconsistent(ac::AffineConstraints,idx,dim)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac 
    return any( hi(bound) < 0 && last_non_zero < idx for (bound,last_non_zero) in zip(c_values,c_last_non_zero_coordinates))   
end

function update_constraints(ac::AffineConstraints,idx::Int,val_at_idx::boxed_nf_elem)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac
    
    return (c_vectors,[b - val_at_idx*r[idx] for (r,b) in zip(c_vectors,c_values)],c_last_non_zero_coordinates)

end


function interval_for_k_j(field,ac::AffineConstraints,j::Int)::Interval{Float64}
    applicable = [(vec[j],val) for (vec,val,last) in zip(ac...) if last == j]
    pos = [v//c for (c,v) in applicable if lo(c) > 0]
    neg = [v//c for (c,v) in applicable if hi(c) < 0]
    no_lb = isempty(neg)
    no_ub = isempty(pos)
    lb = ( no_lb ? -∞  : maximum(x->lo(x),neg) )::Float64
    ub = ( no_ub ? ∞ : minimum(x->hi(x),pos) )::Float64
    return IntervalArithmetic.is_valid_interval(lb,ub) ? lb..ub : ∅
end

# this looks like it should be faster (no temporary arrays created) but it results in slower runtime…
function interval_for_k_j2(field,ac::AffineConstraints,j::Int)::Interval{Float64}
    lb = -∞::Float64
    ub = ∞::Float64
    for (vec,val,last) in zip(ac...) 
        if last == j
            if lo(vec[j])>0
                ub = max(ub,lo(val//vec[j]))
            end
            if hi(vec[j]) < 0
                lb = min(lb,hi(val//vec[j]))
            end
        end
    end
    return IntervalArithmetic.is_valid_interval(lb,ub) ? lb..ub : ∅
end



@inline function _add_if_valid_root(
    vd::VinbergData,
    stem::Vector{boxed_nf_elem},
    stem_can_rep::Vector{boxed_nf_elem},
    stem_length::Int,
    root_length_idx::Int,
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

    stem_no_box = ex.(stem)
    stem_can_no_box = ex.(stem_can_rep)

    @toggled_assert stem_can_no_box == to_can_rep(vd,stem_no_box) "Sanity check. `stem_can_rep` must always be equal to `to_can_rep(vd,stem)`!"
    @toggled_assert times(vd,stem_can_no_box,stem_can_no_box) == l "Sanity check. If we are here, by the previous case (j==vd.dim) necessarily we have the right length."
    if !all(ex(bound) ≥ 0 for bound in c_values)
        return 
    end

    if is_integral(space, ring, stem_can_no_box) && is_root(space,ring,stem_can_no_box,l) 
        # integralness should be guaranteed to hold for all but the last coordinate I think.
        append!(roots,[deepcopy(stem_can_no_box)])
    end
end


# copied from julia Base along with the two functions below, then modified
midpoint(lo::Int, hi::Int) = lo + ((hi - lo) >>> 0x01)

# Ideally for each element sk we'd store both upper and lower bound and could compare x_lower to sk_upper here in a cheap but "secure" fashion (our bounds might be a bit loose, but that wouldn't matter).
# Sadly when I try to store both upper and lower in the T2Cache, something somewhere kills performance, and I could not find what. This is in 7f3d4d5b5f240568aa536114d578f728f0095152 and should be researched further!
function custom_searchsortedfirst(v, x)
    u = 1
    lo = 0 
    hi = length(v)+1
    @inbounds while lo < hi - u
        m = midpoint(lo, hi)
        if ap(v[m][1]).hi < x
            lo = m
        else
            hi = m
        end
    end
    return hi
end

# index of the last value of vector a that is less than or equal to x;
# returns 0 if x is less than all values of v.
function custom_searchsortedlast(v, x)
    u = 1
    lo = 0
    hi = length(v)+1
    @inbounds while lo < hi - u
        m = midpoint(lo, hi)
        if x < ap(v[m][1]).lo
            hi = m
        else
            lo = m
        end
    end
    return lo
end

@inline function find_range(
    field,
    interval_αk::Interval,
    α_over_s::boxed_nf_elem,
    ordered::Vector
   )::Tuple{Int,Int,Int,Int}
   
    lb  = (interval_αk/ap(α_over_s)).lo 
    ub  = (interval_αk/ap(α_over_s)).hi
    no_lb = lb == -∞
    no_ub = ub == ∞
    
    if lb ≤ 0 && ub ≥ 0 

        last_idx_neg = (no_lb ? length(ordered) : custom_searchsortedlast(ordered,-lb) )
        last_idx_pos = (no_ub ? length(ordered) : custom_searchsortedlast(ordered,ub) )
        @toggled_assert all(ap(x[1]).lo > ub for x in ordered[last_idx_pos+1:end])
        @toggled_assert all(ap(-x[1]).hi < lb for x in ordered[last_idx_neg+1:end])
        
        return (1,last_idx_pos,1,last_idx_neg)

    elseif lb ≥ 0 && ub ≥ 0
        
        first_idx_pos = (no_lb ? 1 : custom_searchsortedfirst(ordered,lb) )
        last_idx_pos = (no_ub ? length(ordered) : custom_searchsortedlast(ordered,ub) )
        @toggled_assert all(ap(x[1]).lo >ub for x in ordered[last_idx_pos+1:end])
        @toggled_assert all(ap(x[1]).hi <lb for x in ordered[1:first_idx_pos-1])

        return (first_idx_pos,last_idx_pos,1,0)


    elseif lb ≤ 0 && ub ≤ 0
        
        first_idx_neg = (no_ub ? 1 : custom_searchsortedfirst(ordered,-ub) )
        last_idx_neg = (no_lb ? length(ordered) : custom_searchsortedlast(ordered,-lb) )
        @toggled_assert all(ap(-x[1]).lo >ub for x in ordered[1:first_idx_neg-1])
        @toggled_assert all(ap(-x[1]).hi <lb for x in ordered[last_idx_neg+1:end])
        
        return (1,0,first_idx_neg,last_idx_neg)

    end
    
    @assert false "should not be reachable"

end




const ExactInterval = Tuple{Tuple{Bool,nf_elem},Tuple{Bool,nf_elem}}

function exact_interval_for_k_j(field,ac::AffineConstraints,j::Int)::ExactInterval
    applicable = [(ex(vec[j]),ex(val)) for (vec,val,last) in zip(ac...) if last == j]
    pos = [v//c for (c,v) in applicable if c > 0]
    neg = [v//c for (c,v) in applicable if c < 0]
    no_lb = isempty(neg)
    no_ub = isempty(pos)
    lb = ( no_lb ? field(1) : maximum(neg) )::nf_elem
    ub = ( no_ub ? field(0) : minimum(pos) )::nf_elem
    return ((no_lb,lb),(no_ub,ub))
end

function in_interval(
    k::boxed_nf_elem,
    interval::ExactInterval
)
    ((no_lb,lb),(no_ub,ub)) = interval
    return (no_lb || ex(k) ≥ lb) && (no_ub || ex(k) ≤ ub)
end

function is_empty(i::ExactInterval)
    (no_lb,lb),(no_ub,ub) = i
    return !no_lb && !no_ub && (lb > ub)
end

function exact_interval_in_approx_interval(ei::ExactInterval,ai::Interval)
    ((no_lb,lb),(no_ub,ub)) = ei
    
    #@info "exact interval = $ei ~ [$(lower_float64(lb)),$(upper_float64(ub))]"
    #@info "approx interval = $ai"

    no_lb && ai.lo ≠ -∞ && return false
    no_ub && ai.hi ≠ ∞ && return false
    ai == ∅ && !is_empty(ei) && return false
    is_empty(ei) && return true
    ok =  (upper_float64(lb) ≥ ai.lo) && (lower_float64(ub) ≤ ai.hi)
    if !ok
        @error "exact interval = $ei ~ [$(lower_float64(lb)),$(upper_float64(ub))] is empty? $(is_empty(ei))"
        @error "approx interval = $ai"
        @error "-------------------------------------"
        @error "$(upper_float64(lb))≥ $(ai.lo) ? $(upper_float64(lb) ≥ ai.lo)"
        @error "$(lower_float64(ub))≤ $(ai.hi) ? $(lower_float64(ub) ≥ ai.hi)"
    end
    return ok
end


function _extend_root_stem!(
    vd::VinbergData,
    stem::Vector{boxed_nf_elem},
    stem_can_rep::Vector{boxed_nf_elem},
    stem_length::Int,
    root_length_idx::Int,
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
    
    if root_length_minus_stem_norm_squared == 0 # no more space in the root_length. I don't think this present "short_circuiting" necessarily helps a lot but it might do
        return _add_if_valid_root(vd,stem,stem_can_rep,vd.dim,root_length_idx,root_length,root_length_minus_stem_norm_squared,constraints,t2_cache,roots)
    end

    if j == vd.dim + 1
        return _add_if_valid_root(vd,stem,stem_can_rep,stem_length,root_length_idx,root_length,root_length_minus_stem_norm_squared,constraints,t2_cache,roots)
    end

    #if j == vd.dim
    #    return _extend_root_stem_one_coord_left!(vd,stem,stem_can_rep,stem_length,root_length_idx,root_length,root_length_minus_stem_norm_squared,constraints,t2_cache,roots)
    #end

    field = vd.field
    ring = vd.ring
    space = vd.quad_space
    P = infinite_places(field)
    l = root_length
    @toggled_assert vd.possible_root_norms_squared_up_to_squared_units[root_length_idx] == l
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
    two_α_over_ls = vd.two_diago_over_scaling_times_length[root_length_idx][j]


  
    l_j_ssq_over_α = l_j//ex(α_over_s²)
    
    #t2_bound_for_sk = approx_sum_at_places(l_j//(α_over_s²),first_place_idx=1)+1
    t2_bound_for_sk = Hecke.tr(l_j_ssq_over_α)
    last_bounded_t2_candidates_vector_idx = bounded_t2_elems!(
        vd,
        t2_bound_for_sk,
        t2_cache
    )
  

    #@toggled_assert l_j_ssq_over_α == 0 || all(ispositive(l_j_ssq_over_α,p) for p in P)
    conjugates_bound = sqrt.((get_enclosing_intervals_float64(l_j_ssq_over_α)))
    
    #    Constraint on norm at al places:
    #
    #       k^2α ≤ l_j at all places
    #    ⇔  (sk)^2 α/s^2 ≤ l_j at all places
    #    ⇔  (sk)^2 * α_over_s² ≤ l_j at all places
    #    ⇔  |(sk)| ≤ √(l_j//α_over_s²) at all places
    #
    #
    exact_good_norm(sk) = all( ≤(sk^2,l_j_ssq_over_α,p) for p in P) # use this one if need to be sure
    function good_norm(sk,boxes)
        is_in_bound = all(!(σbound.hi < min(abs(σsk.lo),abs(σsk.hi)))  for (σsk,σbound) in zip(boxes,conjugates_bound)) 
        # σbound ≺ abs(σsk) is would be cleaner (note that the symbol is `\prec`) but waiting for https://github.com/JuliaIntervals/IntervalArithmetic.jl/pull/443 
        if Options.conjugates_bound_check_mode() == Options.OnlyFloaty
            return is_in_bound 

        elseif Options.conjugates_bound_check_mode() == Options.FloatyPlusExactCheck
            @assert !is_in_bound ⇒ !exact_good_norm(ex(sk))
            return is_in_bound
            # Normally, the way we approximated the elements in abs_conjugates_sk and conjugates_bound (lower and upper respectively), we should be safe
            # That is, if `is_in_bound` doesn't hold, then we really are out of bounds, so no good root lost.
        elseif Options.conjugates_bound_check_mode() == Options.Exact
            return exact_good_norm(sk)
        end
    end

    #    Constraint as given by the crystallographic condition:
    #
    #       l | 2kα  
    #    ⇔  2kα/l ∈ ring  
    #    ⇔  sk (2α//ls) ∈ ring 
    #    ⇔  sk (two_α_over_ls) ∈ ring 
    crystal(sk,mat) = begin
        @toggled_assert (sk * ex(two_α_over_ls) ∈ ring)  ⇔  mat[root_length_idx,j]
        return mat[root_length_idx,j] 
    end

    integral(a_stem_can_rep) = all((ex(a_stem_can_rep[idx]) ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j])

    
    # The idea is that interval_k_j gives an interval outside of which k_jα_j is not valid due to the constraints of acute angles given by previous roots.
    # The code below SHOULD then use this interval to only iterate over k_js in this interval.
    interval_αk = interval_for_k_j(field,constraints,j) 
    #interval_αk = -∞..∞
    #exact_interval_αk = exact_interval_for_k_j(field,constraints,j)
    #@assert exact_interval_in_approx_interval(exact_interval_αk,interval_αk)

    # If the endpoints are not ±∞, rescale them to get endpoints for sk instead of endpoints of α*k
    @assert ex(α_j) > 0
    if interval_αk == ∅
        return 
    end

    stem_updated = deepcopy(stem)
    stem_can_rep_updated = deepcopy(stem_can_rep) #copy(stem_can_rep)

    last_coordinate = j == vd.dim
    for (idx,ordered) in enumerate(last_coordinate ? (t2_cache.elems[last_bounded_t2_candidates_vector_idx:last_bounded_t2_candidates_vector_idx]) : (t2_cache.elems[1:last_bounded_t2_candidates_vector_idx]) )
       
        @toggled_assert issorted(ordered,by=x->ex(x[1]))
        isempty(ordered) && continue
        
        check_t2 = idx == last_bounded_t2_candidates_vector_idx
     
        (first_idx_pos,last_idx_pos,first_idx_neg,last_idx_neg) = find_range(field,interval_αk, α_over_s, ordered) 
        for i in min(first_idx_pos,first_idx_neg):max(last_idx_pos,last_idx_neg)
            
            sk,boxes,t2sk,crystal_mat = ordered[i]
            pos =  first_idx_pos ≤ i && last_idx_pos ≥ i
            neg = first_idx_neg ≤ i && last_idx_neg ≥ i

            if ( !last_coordinate ⇒ ( (check_t2 ⇒ (t2sk ≤ t2_bound_for_sk)) && good_norm(sk,boxes) && crystal(ex(sk),crystal_mat) ) ) &&
                ( last_coordinate ⇒ (ex(sk)^2 == l_j_ssq_over_α && crystal(ex(sk),crystal_mat) ))
                #sk * two_α_over_ls ∈ ring && # crystallographic condition 
                #all( ≤(sk^2 * α_over_s²,l_j,p) for p in P) #  norms are OK
                
                k = sk // s_j
                if pos 
                    @toggled_assert !(ap(k*α_j) ≻ interval_αk) && !(ap(k*α_j) ≺ interval_αk)
                    stem_updated = copy!(stem_updated,stem); stem_updated[j] = k
                    #stem_can_rep_updated = [stem_can_rep[i] + k * box(v_j[i]) for i in 1:vd.dim]
                    u_plus_k_v(stem_can_rep_updated,stem_can_rep,k,v_j)
                    if integral(stem_can_rep_updated) # all((stem_can_rep_updated[idx] ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j]) # integral
                        _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,root_length_idx,l,l_j - ex(k)^2*ex(α_j),update_constraints(constraints,j,k*α_j),t2_cache,roots)
                    end
                end
                if neg && ex(k)≠0 
                    @toggled_assert !(ap(-k*α_j) ≻ interval_αk) && !(ap(-k*α_j) ≺ interval_αk)
                    stem_updated = copy!(stem_updated,stem); stem_updated[j] = -k
                    #stem_can_rep_updated = [stem_can_rep[i] - k * box(v_j[i]) for i in 1:vd.dim]
                    u_plus_k_v(stem_can_rep_updated,stem_can_rep,-k,v_j)
                    if integral(stem_can_rep_updated) # all((stem_can_rep_updated[idx] ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j]) # integral
                        _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,root_length_idx,l,l_j - ex(k)^2*ex(α_j),update_constraints(constraints,j,-k*α_j),t2_cache,roots)
                    end
                end

            end

        end
    end



end

≻(i1,i2) = i2 ≺ i1
function u_plus_k_v(to,u,k,v)
    @assert length(to) == length(u)
    @assert length(to) == length(v)
    @inbounds for i in 1:length(to)
        mul!(to[i].x,k.x,v[i].x)
        addeq!(to[i].x,u[i].x)
        to[i].box = u[i].box + k.box * v[i].box
    end
    #@assert to = u .+ (k .* v)
end

function extend_root_stem(
    vd::VinbergData,
    stem_diag_rep::Vector{nf_elem},
    root_length,
    constraints::AffineConstraints,
    t2_cache::BoundedT2ElemsCache
)
    
    root_length_idx = findfirst(x->x==root_length,vd.possible_root_norms_squared_up_to_squared_units)
    
    stem_length = length(stem_diag_rep)
    stem_diag_rep = (vcat(stem_diag_rep,fill(vd.field(0),vd.dim-length(stem_diag_rep))))
    stem_can_rep = (to_can_rep(vd,stem_diag_rep))
    
    stem_norm_squared = norm_squared(vd,stem_can_rep)
    
    roots_go_here = Vector{Vector{nf_elem}}()
    _extend_root_stem!(vd,box.(stem_diag_rep),box.(stem_can_rep),stem_length,root_length_idx,ex(root_length),ex(root_length)-stem_norm_squared,constraints,t2_cache,roots_go_here)
     return roots_go_here
end

function roots_at_distance_zero(vd::VinbergData)

    t2_cache = BoundedT2ElemsCache(vd) 
    zero_stem = nf_elem[vd.field(0)]
    
    return vcat([extend_root_stem(vd,zero_stem,l,no_constraints(),t2_cache) for l in vd.possible_root_norms_squared_up_to_squared_units]...)
end

function cone_roots(vd,roots_at_distance_zero)

    @warn "Cone roots computation are approximative ⇒ double check the results by hand."
    @warn "If results look dubious, increase LP precision by calling `VinbergsAlgorithm.Options.set_lp_precision(desired_precision::Int)` (currently $(Options.lp_precision()))"
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
        if  all((-1)*r ≠ cr for cr in cone_roots)
            if is_necessary_halfspace(vd.gram_matrix.entries,cone_roots,-r)
                push!(cone_roots,r)
            end
        
        end
    end
   
    cone_roots = drop_redundant_halfspaces(vd.gram_matrix.entries,cone_roots)
    @assert all(times(vd,r₁,r₂) ≤ 0 for r₁ in cone_roots for r₂ in cone_roots if r₁≠r₂)

    for r in roots_at_distance_zero
        @toggled_assert !is_necessary_halfspace(vd.gram_matrix.entries, cone_roots,r) || !is_necessary_halfspace(vd.gram_matrix.entries, cone_roots,-r)
    end

    for r in cone_roots
        @info r
    end
    return cone_roots

end

function cone_roots(vd)
    cone_roots(vd,roots_at_distance_zero(vd))
end

function roots_for_pair(vd,pair,prev_roots;t2_cache=nothing)
    
    @info "roots_for_pair(vd,$pair,$(length(prev_roots)))"

    if t2_cache === nothing
        t2_cache = BoundedT2ElemsCache(vd) 
    end

    (k,l) = pair
    fake_dist = -(k^2)*ex(vd.diagonal_values[1])//ex(l)
   
    ######################## The following only useful for non-diags
    all_zero_on_coord_after(vd,vector_idx,coord_idx) = all(vd.diagonal_basis[l][coord_idx] == 0 for l in vector_idx+1:vd.dim)
    if !all(all_zero_on_coord_after(vd,1,idx) ⇒ (k*ex((vd.diagonal_basis[1][idx])) ∈ vd.ring) for idx in 1:vd.dim)
        @info "$k would define non-integral coordinates ⇒ dropping it"
        return [] 
    end
    ######################
    
    # Construct the constraints defined by the previous roots
    prev_roots_diag = [to_diag_rep(vd,prev_root) for prev_root in prev_roots]
    prev_roots_diag_bound = [-k*ex(vd.diagonal_values[1])*prev_root[1] for prev_root in prev_roots_diag]
    prev_roots_constraints = MkAffineConstraints([box.(r) for r in prev_roots_diag],box.(prev_roots_diag_bound))

    #@info "roots_for_pair($pair,$prev_roots)"
    roots = extend_root_stem(vd,[k],l,prev_roots_constraints,t2_cache)
    
    @assert all(times(vd,root,prev) ≤ 0 for root in roots for prev in prev_roots)  "All angles with previous roots should be acute."
    @assert all(is_root(vd,root) for root in roots) "All outputs of extend_root_stem must be roots"
    @assert all(norm_squared(vd,root) == ex(l) for root in roots) "All outputs of extend_root_stem must have correct length"
    @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
    @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂) "Two roots at same distance and compatible with everything before them should be compatible with each other." 
   
    
    return roots
    
end

function roots_for_next_pair!(vd,dict,prev_roots;t2_cache=nothing)

    pair = next_min_pair!(vd,dict)
    k,l = pair
    fake_dist = -(k^2)*ex(vd.diagonal_values[1])//ex(l)
    @info "next pair is $(k),$(ex(l)) (fake_dist is $fake_dist ≈ $(approx(fake_dist,8)))"


    # Here we should assert that the roots in prev_roots are actually closer than the next min pair can give us, otherwise we will get inconsistencies

    roots =  roots_for_pair(vd,pair,prev_roots,t2_cache=t2_cache)

    @toggled_assert all(fake_dist_to_basepoint(vd,r) == fake_dist for r in roots)
    return roots

end

function roots_for_next_pairs!(vd,dict,prev_roots;t2_cache=nothing)
   
    nthreads = Threads.nthreads()
    #nthreads = 1
    pairs = [next_min_pair!(vd,dict) for _ in 1:nthreads]
    rootss = [Vector{nf_elem}[] for _ in 1:nthreads]
    Threads.@threads for (i,pair) in collect(enumerate(pairs))
        (k,l) = pair
        fake_dist = -(k^2)*ex(vd.diagonal_values[1])//ex(l)
        @info "[$i]next pair is $(k),$(ex(l)) (fake_dist is $fake_dist ≈ $(approx(fake_dist,8)))"
        rootss[i] = roots_for_pair(vd,pair,prev_roots,t2_cache=t2_cache)
    end
    
    cumulative_roots = rootss[1]
    for i in 2:nthreads
        filter!(r₁ -> all(times(vd,r₁,r₂)≤0 for r₂ in cumulative_roots),rootss[i])
        append!(cumulative_roots,rootss[i])
    end
    return cumulative_roots
end

function next_n_roots!(vd,prev_roots,dict,das;n=10,t2_cache=nothing)

    if t2_cache===nothing
        t2_cache = BoundedT2ElemsCache(vd)
    end
    
    # TODO
    # asserts that the previous roots are closer than what the dict proposes, otherwise we'll get inconsistencies


    roots = prev_roots
    #Coxeter_matrix = get_Coxeter_matrix(vd.quad_space, vd.ring, prev_roots) 
    new_roots = Vector{Vector{nf_elem}}()
    while n > 0 


        new_roots = roots_for_next_pairs!(vd,dict,roots;t2_cache=t2_cache)
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
    das = DiagramAndSubs(Cox_matrix,vd.dim-1)
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
