
struct VinbergData

    dim::Int
    gram_matrix::Matrix{fmpz}

    diagonal_basis::Vector{Vector{fmpz}}
    diagonal_values::Vector{fmpz}
    scaling::Vector{fmpz}
    
    diagonal_change::Array{fmpq,2}
    diagonal_change_inv::Array{fmpq,2}


    possible_root_norms_squared_up_to_squared_units::Vector{fmpz}

    # Precomputed stuff:
    diago_over_scaling::Vector{fmpq}
    diago_over_scalingsq::Vector{fmpq}
    two_diago_over_scaling_times_length::Dict{Int,Vector{fmpq}}
    diago_vector_last_on_coordinates::Vector{Vector{Int}}

end


function VinbergData(gram_matrix,v₀=nothing)

    (n,m) = size(gram_matrix)
    @assert n == m "The Gram gram_matrix must be square."

    #quad_lattice = Hecke.lattice(gram_matrix)
        
    if v₀≠nothing
        @assert v₀' * gram_matrix * v₀ < 0
    end

    diagonal_basis_vecs,diagonal_values,scaling = diagonalize_and_get_scaling(gram_matrix,v₀)
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

    rG = matrix(QQ,gram_matrix)
    display(rG)
    cofactorsG = ZZ.(collect(Hecke.det(rG) * Hecke.inv(rG))) # 
    last_invariant_factor = abs(ZZ(Hecke.det(rG)//gcd(cofactorsG)))
    
    # The possible lengths of roots are the divisor of `2*last_invariant_factor`.
    twice_LIF = 2*last_invariant_factor
    lengths =  [ZZ(k) for k in 1:twice_LIF if twice_LIF%k == 0]
 
    # Precomputations ---
    diago_over_scaling = [α//s for (α,s) in zip(diagonal_values,scaling)]::Vector{fmpq} 
    diago_over_scalingsq = [α//(s^2) for (α,s) in zip(diagonal_values,scaling)]::Vector{fmpq} 
    two_diago_over_scaling_times_length = Dict([idx => (2//l) .* diago_over_scaling for (idx,l) in enumerate(lengths)])::Dict{Int,Vector{fmpq}}
    
    is_last_on_coord(vecs,v_i,c_i) = vecs[v_i][c_i]≠0 && all(vecs[i][c_i] == 0 for i in v_i+1:length(vecs))
    diago_vector_last_on_coordinates =[filter(c_i -> is_last_on_coord(diagonal_basis_vecs,v_i,c_i), 1:n) for v_i in 1:n] ::Vector{Vector{Int}}
    # -------------------

    vd =  VinbergData(
        n,
        gram_matrix,
        [[v for v in vec] for vec in diagonal_basis_vecs],
        [v for v in diagonal_values],
        [v for v in scaling],
        Matrix(matrix(QQ,diagonal_basis_vecs))',    # notice the transpose here and below. Much frustration took place before I found out I needed those!
        Matrix(Hecke.inv(matrix(QQ,diagonal_basis_vecs)))',
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
    @assert (vd.diagonal_values[1]) < 0 "sanity check"

    return (vd.diagonal_basis[1])
end

times(G,u,v) = u' * G * v
times(vd::VinbergData,u,v) = times(vd.gram_matrix,u,v)

norm_squared(matrix,u) = times(matrix,u,u)
norm_squared(vd::VinbergData,u) = times(vd.gram_matrix,u,u)

is_root(vd::VinbergData, vector::Vector) = is_root(vd.gram_matrix,vector)
Gram_coeff(vd::VinbergData,r₁::Vector,r₂::Vector) = Gram_coeff(vd.gram_matrix,r₁,r₂)
Coxeter_coeff(vd::VinbergData,r₁::Vector,r₂::Vector) = Coxeter_coeff(vd.gram_matrix,r₁,r₂)
Coxeter_matrix(vd::VinbergData,roots) = Coxeter_matrix(vd.gram_matrix,roots)

"""
    fake_dist_to_point(vd,point,root)

If `point` defines an element of ``ℍ^n`` and `root` a hyperplane, compute the fake distance between the point and the hyperplane, satisfying
```math
    \\sinh²(dist(point,hyperplane)) = fake_dist(point,hyperplane)
```
"""
function fake_dist_to_point(vd,point,root)
    @toggled_assert is_root(vd.gram_matrix,root) "Need a root"
    @toggled_assert norm_squared(vd,point) < 0 "Need a point in hyperbolic space."

    fake_dist = -times(vd,root,point)^2//(norm_squared(vd,point)*norm_squared(vd,root))

    return fake_dist 
end

function fake_dist_to_basepoint(vd,u)

    

    fake_dist =  fake_dist_to_point(vd,basepoint(vd),u)
    @toggled_assert fake_dist == -(to_diag_rep(vd,u)[1]^2)*(vd.diagonal_values[1])//norm_squared(vd,u) 

    return fake_dist

end





const LeastKByRootNormSquared = Dict{
    fmpq, # norm_squared
    fmpq, #k
}



function next_min_k_for_l(vd,current_k,l)

    s_1 = vd.scaling[1]
    α_1 = vd.diagonal_values[1] 
    
    k = current_k + 1//s_1

    while !divides(l,2*k*α_1)
        k += 1//s_1
    end
    return k
end

function next_min_k_for_l!(vd,dict,l)

    current_k = dict[l] 
    new_k = next_min_k_for_l(vd,current_k,l)
    dict[l] = new_k

end


function init_least_k_by_root_norm_squared(vd::VinbergData)

    return Dict(
        [
            l => next_min_k_for_l(vd,0,l)
            for l in vd.possible_root_norms_squared_up_to_squared_units
        ]
    )
end


function next_min_pair(vd,dict)
    min_pair = nothing

    # that's the value we want to minimize (related to the distance)
    val(k,l) = k^2//l
    

    for (l,k) in dict 
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


const AffineConstraints = Tuple{
    Vector{Vector{fmpq}}, # Previous roots
    Vector{fmpq},         # ≤ bound
    Vector{Int},           # last non non-zero coordinate
}     


function MkAffineConstraints(
    vecs::Vector{Vector{fmpq}}, # Previous roots
    vals::Vector{fmpq},         # ≤ bound
)
    @assert all(any(x->x≠0,vec) for vec in vecs)
    last_nz = Int[findlast(x->x≠0,vec) for vec in vecs]
    return (vecs,vals,last_nz)
end

no_constraints() = MkAffineConstraints(Vector{Vector{fmpq}}(),Vector{fmpq}())

function clearly_inconsistent(ac::AffineConstraints,idx,dim)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac 
    return any( bound < 0 && last_non_zero < idx for (bound,last_non_zero) in zip(c_values,c_last_non_zero_coordinates))   
end

function update_constraints(ac::AffineConstraints,idx::Int,val_at_idx::fmpq)
    (c_vectors,c_values,c_last_non_zero_coordinates) = ac
    
    return (c_vectors,[b - val_at_idx*r[idx] for (r,b) in zip(c_vectors,c_values)],c_last_non_zero_coordinates)

end




@inline function _add_if_valid_root(
    vd::VinbergData,
    stem::Vector{fmpq},
    stem_can_rep::Vector{fmpq},
    stem_length::Int,
    root_length_idx::Int,
    root_length::fmpq,
    root_length_minus_stem_norm_squared::fmpq,
    constraints::AffineConstraints,
    roots::Vector{Vector{fmpq}}
)

    space = vd.gram_matrix
    l = root_length   
    (c_vectors,c_values,c_last_non_zero_coordinates) = constraints


    @toggled_assert stem_can_rep == to_can_rep(vd,stem) "Sanity check. `stem_can_rep` must always be equal to `to_can_rep(vd,stem)`!"
    
    @toggled_assert times(vd,stem_can_rep,stem_can_rep) == l

    if !all(bound ≥ 0 for bound in c_values)
        return 
    end

    if is_integral(stem_can_rep) && is_root(space,stem_can_rep,l) 
        # integralness should be guaranteed to hold for all but the last coordinate I think.
        append!(roots,[deepcopy(stem_can_rep)])
    end
end


const ExactInterval = Tuple{Tuple{Bool,fmpq},Tuple{Bool,fmpq}}

function exact_interval_for_k_j(ac::AffineConstraints,j::Int)::ExactInterval
    applicable = [((vec[j]),(val)) for (vec,val,last) in zip(ac...) if last == j]
    pos = [v//c for (c,v) in applicable if c > 0]
    neg = [v//c for (c,v) in applicable if c < 0]
    no_lb = isempty(neg)
    no_ub = isempty(pos)
    lb = ( no_lb ? QQ(1) : maximum(neg) )::fmpq
    ub = ( no_ub ? QQ(0) : minimum(pos) )::fmpq
    return ((no_lb,lb),(no_ub,ub))
end

function in_interval(
    k::fmpq,
    interval::ExactInterval
)
    ((no_lb,lb),(no_ub,ub)) = interval
    return (no_lb || (k) ≥ lb) && (no_ub || (k) ≤ ub)
end

function is_empty(i::ExactInterval)
    (no_lb,lb),(no_ub,ub) = i
    return !no_lb && !no_ub && (lb > ub)
end




function _extend_root_stem!(
    vd::VinbergData,
    stem::Vector{fmpq},
    stem_can_rep::Vector{fmpq},
    stem_length::Int,
    root_length_idx::Int,
    root_length::fmpq,
    root_length_minus_stem_norm_squared::fmpq,
    constraints::AffineConstraints,
    roots::Vector{Vector{fmpq}}
)
    
    j = stem_length + 1 
    
    if clearly_inconsistent(constraints,j,vd.dim)
        println("clearly inconsistent")
        return
    end
    
    if root_length_minus_stem_norm_squared == 0 # no more space in the root_length. I don't think this present "short_circuiting" necessarily helps a lot but it might do
        return _add_if_valid_root(vd,stem,stem_can_rep,vd.dim,root_length_idx,root_length,root_length_minus_stem_norm_squared,constraints,roots)
    end

    if j == vd.dim + 1
        return _add_if_valid_root(vd,stem,stem_can_rep,stem_length,root_length_idx,root_length,root_length_minus_stem_norm_squared,constraints,roots)
    end

    #if j == vd.dim
    #    return _extend_root_stem_one_coord_left!(vd,stem,stem_can_rep,stem_length,root_length_idx,root_length,root_length_minus_stem_norm_squared,constraints,roots)
    #end

    space = vd.gram_matrix
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


  
    l_j_ssq_over_α = l_j//α_over_s²
    

  


    #    Constraint as given by the crystallographic condition:
    #
    #       l | 2kα  
    #    ⇔  2kα/l ∈ ring  
    #    ⇔  sk (2α//ls) ∈ ring 
    #    ⇔  sk (two_α_over_ls) ∈ ring
    crystal(k) = divides(l,2*k*α_j)

    integral(a_stem_can_rep) = all(isinteger(a_stem_can_rep[idx]) for idx in vd.diago_vector_last_on_coordinates[j])

    
    interval_αk = exact_interval_for_k_j(constraints,j) 
    if is_empty(interval_αk) == ∅
        println("empty interval")
        return 
    end

    stem_updated = deepcopy(stem)
    stem_can_rep_updated = deepcopy(stem_can_rep) #copy(stem_can_rep)

    
    step1 = l//(2α_j)
    step2 = 1//s_j
    C = denominator(step1) * denominator(step2)
    step = lcm(ZZ(step1*C),ZZ(step2*C))//C

    k = fmpq(0)
    while k^2 ≤ l_j//α_j
        
        #= 
        if j == vd.dim
            println("-------------------------------")
            println("l = $l")
            println("l_j = $l_j")
            println()
        end
        =#
        #@assert crystal(k)
        if ((j == vd.dim) ⇒ (k^2 == l_j//α_j)) 
            
            kα_j = k*α_j
            if in_interval(kα_j,interval_αk)
                stem_updated = copy!(stem_updated,stem); stem_updated[j] = k
                stem_can_rep_updated = [stem_can_rep[i] + k * (v_j[i]) for i in 1:vd.dim]
                if integral(stem_can_rep_updated) # all((stem_can_rep_updated[idx] ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j]) # integral
                    _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,root_length_idx,l,l_j - k^2*(α_j),update_constraints(constraints,j,kα_j),roots)
                end
            end
            if in_interval(-kα_j,interval_αk) && (k)≠0 
                stem_updated = copy!(stem_updated,stem); stem_updated[j] = -k
                stem_can_rep_updated = [stem_can_rep[i] - k * (v_j[i]) for i in 1:vd.dim]
                if integral(stem_can_rep_updated) # all((stem_can_rep_updated[idx] ∈ ring) for idx in vd.diago_vector_last_on_coordinates[j]) # integral
                    _extend_root_stem!(vd,stem_updated,stem_can_rep_updated,j,root_length_idx,l,l_j - (k)^2*(α_j),update_constraints(constraints,j,-kα_j),roots)
                end
            end

        else
            # nothing to see here
        end
        
        k += step #1//s_j
    end
    


end

function extend_root_stem(
    vd::VinbergData,
    stem_diag_rep::Vector{fmpq},
    root_length,
    constraints::AffineConstraints,
)
    
    root_length_idx = findfirst(x->x==root_length,vd.possible_root_norms_squared_up_to_squared_units)
    
    stem_length = length(stem_diag_rep)
    stem_diag_rep = (vcat(stem_diag_rep,fill(QQ(0),vd.dim-length(stem_diag_rep))))
    stem_can_rep = (to_can_rep(vd,stem_diag_rep))
    
    stem_norm_squared = norm_squared(vd,stem_can_rep)
    
    roots_go_here = Vector{Vector{fmpq}}()
    _extend_root_stem!(vd,(stem_diag_rep),(stem_can_rep),stem_length,root_length_idx,QQ(root_length),(root_length)-stem_norm_squared,constraints,roots_go_here)
     return roots_go_here
end

function roots_at_distance_zero(vd::VinbergData)

    zero_stem = fmpq[0]
    
    return vcat([extend_root_stem(vd,zero_stem,l,no_constraints()) for l in vd.possible_root_norms_squared_up_to_squared_units]...)
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
    
    cone_roots = Vector{Vector{fmpq}}()


    for r in roots_at_distance_zero
        @debug "looking at $r"
        if  all((-1)*r ≠ cr for cr in cone_roots)
            @debug "so far so good"
            if is_necessary_halfspace(vd.gram_matrix,cone_roots,-r)
                @debug "degeneration"
                push!(cone_roots,r)
            end
        
        end
        @debug "have $(length(cone_roots)) cone roots" 
    end
    
    cone_roots = drop_redundant_halfspaces(vd.gram_matrix,cone_roots)
    @assert all(times(vd,r₁,r₂) ≤ 0 for r₁ in cone_roots for r₂ in cone_roots if r₁≠r₂)

    for r in cone_roots
        @info r
    end
    return cone_roots

end

function cone_roots(vd)
    cone_roots(vd,roots_at_distance_zero(vd))
end

function roots_for_pair(vd,pair,prev_roots)
    


    (k,l) = pair
    fake_dist = -(k^2)*(vd.diagonal_values[1])//(l)
   
    ######################## The following only useful for non-diags
    all_zero_on_coord_after(vd,vector_idx,coord_idx) = all(vd.diagonal_basis[l][coord_idx] == 0 for l in vector_idx+1:vd.dim)
    if !all(all_zero_on_coord_after(vd,1,idx) ⇒ isinteger(k*((vd.diagonal_basis[1][idx]))) for idx in 1:vd.dim)
        @info "$k would define non-integral coordinates ⇒ dropping it"
        return [] 
    end
    ######################
    
    # Construct the constraints defined by the previous roots
    prev_roots_diag = [to_diag_rep(vd,prev_root) for prev_root in prev_roots]
    prev_roots_diag_bound = [-k*(vd.diagonal_values[1])*prev_root[1] for prev_root in prev_roots_diag]
    prev_roots_constraints = MkAffineConstraints([(r) for r in prev_roots_diag],(prev_roots_diag_bound))

    #@info "roots_for_pair($pair,$prev_roots)"
    roots = extend_root_stem(vd,[k],l,prev_roots_constraints)
    
    @assert all(times(vd,root,prev) ≤ 0 for root in roots for prev in prev_roots)  "All angles with previous roots should be acute."
    @assert all(is_root(vd,root) for root in roots) "All outputs of extend_root_stem must be roots"
    @assert all(norm_squared(vd,root) == (l) for root in roots) "All outputs of extend_root_stem must have correct length"
    @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
    @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂) "Two roots at same distance and compatible with everything before them should be compatible with each other." 
   
    
    return roots
    
end

function roots_for_next_pair!(vd,dict,prev_roots)

    pair = next_min_pair!(vd,dict)
    k,l = pair
    fake_dist = -(k^2)*(vd.diagonal_values[1])//(l)
    @info "next pair is $(k),$((l)) (fake_dist is $fake_dist ≈ $(fake_dist))"


    # Here we should assert that the roots in prev_roots are actually closer than the next min pair can give us, otherwise we will get inconsistencies

    roots =  roots_for_pair(vd,pair,prev_roots)

    @toggled_assert all(fake_dist_to_basepoint(vd,r) == fake_dist for r in roots)
    return roots

end

function next_n_roots!(vd,prev_roots,dict,das;n=10)


    
    # TODO
    # asserts that the previous roots are closer than what the dict proposes, otherwise we'll get inconsistencies


    roots = prev_roots
    #Coxeter_matrix = get_Coxeter_matrix(vd.gram_matrix, vd.ring, prev_roots) 
    new_roots = Vector{Vector{fmpq}}()
    while n > 0 


        new_roots = roots_for_next_pair!(vd,dict,roots)
        n = n - length(new_roots)

        for root in new_roots 
            @info "Got new root $root"
            extend!(das,[Coxeter_coeff(vd.gram_matrix, old,root) for old in roots])
            push!(roots,root)
        end


        if CoxeterDiagrams.is_finite_volume(das)
            # Sanity checks
            @assert all(is_root(vd.gram_matrix,r) for r in roots)
            @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂)
            @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
            @assert all(fake_dist_to_basepoint(vd,roots[i]) ≤ fake_dist_to_basepoint(vd,roots[i+1]) for i in 1:length(roots)-1)
            return (true,(roots,dict,das))
        end
    end

    # Sanity checks
    @assert all(is_root(vd.gram_matrix,r) for r in roots)
    @assert all(times(vd,r₁,r₂)≤0 for r₁ in roots for r₂ in roots if r₁≠r₂)
    @assert all(times(vd,r,basepoint(vd)) ≤ 0 for r in roots) "All outputs must have the basepoint on their negative side."
    @assert all(fake_dist_to_basepoint(vd,roots[i]) ≤ fake_dist_to_basepoint(vd,roots[i+1]) for i in 1:length(roots)-1)

    return (false,(roots,dict,das))
end

function next_n_roots!(
    vd::VinbergData,
    prev_roots::Vector;
    n=10,
)

    # TODO: make this work even when the previous roots are not cone roots.
    # In this case, one just need to iterate over the pairs (k,l) in the dictionary until we're at distance ≥ than the farthest root in prev_roots
    @assert all(fake_dist_to_basepoint(vd,r) == 0 for r in prev_roots) "Only works with previous roots == cone roots"

    Cox_matrix = Coxeter_matrix(vd.gram_matrix, prev_roots) 
    das = CoxeterDiagrams.DiagramAndSubs(Cox_matrix,vd.dim-1)
    dict = init_least_k_by_root_norm_squared(vd)
    
    return next_n_roots!(vd,prev_roots,dict,das;n=n)
end

function next_n_roots!(
    vd::VinbergData;
    n=10,
)
    roots = cone_roots(vd)

    return next_n_roots!(vd,roots;n=n)
end

function all_in_one(
    diagonal::Vector{fmpq},
    n::Int
)
    @assert n>0 && n≤100 "Need a reasonable number of roots to look for."
    @assert length(diagonal) > 2 "Need to have a long enough diagonal."
    @assert all(parent(v) == parent(diagonal[1]) for v in diagonal) "The elements of the diagonal must lie in a common field."
    

    return next_n_roots!(VinbergData(diagm(diagonal)),n=n)
end
