# TODO
#
#=
# 
# candidate is a tuple (indices,gram_matrix,lightgraph,vector_matrix,inverse_of_vector_matrix)
#
#
#
# =#


function is_infinite_order_isometry(
    vd,
    t,
    check_invertible=false,
    check_integral=false,
    check_integrally_invertible=false,
    check_preserves_form=false,
    check_preserves_upper_sheet=false,
)

    # Normally: if it preserve the form and is invertible and integral, the inverse also preserves the form
    # Also, having unit determinant is enough to have an integral inverse
    # But we check too much: better safe than sorry

    check_invertible && @assert isinvertible(t)
    check_integral && @assert is_integral(vd,t)
    check_integrally_invertible && @assert isunit(vd.ring(det(t)))
    check_integrally_invertible && @assert is_integrally_invertible(vd,t)
    check_preserves_form && @assert preserves_the_form(vd,t) 
    check_preserves_upper_sheet && @assert preserves_upper_sheet(vd,t) 

   rk_fixed_t,fixed_t = nullspace(t - identity_matrix(vd.field,vd.dim))

    if rk_fixed_t ≠ 0

        gram_fixed_t = [c' * vd.gram_matrix.entries * d for c in eachcol(fixed_t.entries), d in eachcol(fixed_t.entries)]
        
        d,p = diagonalize_in_field(vd.field,gram_fixed_t)
        
        if any(v < 0 for v in LinearAlgebra.diag(d))
            return false 
        else
            return true
        end

    else
        return true
    end


end


function inf_ord_sym(vd,roots,das,type=:n_walls)
    if type === :n_plus_one_walls
        return inf_ord_sym_n_plus_one_walls(vd,roots,das)
    else
        return inf_ord_sym_n_walls(vd,roots,das)
    end
end

function inf_ord_sym_n_plus_one_walls(vd,roots,das)
    
    Gram = [c' * vd.gram_matrix.entries * d for c in roots, d in roots]
    Gram_hash = hash.(Gram) 


    grouped_candidates = Dict{
                              Vector{Vector{UInt64}},
                              Vector{Tuple{Vector{Int64}, Matrix{nf_elem}, Tuple{LightGraphs.SimpleGraphs.SimpleGraph{Int64}, Dict{Any, Any}, Dict{Any, Any}}}}
                             }()
    for labels in  Combinatorics.powerset(collect(1:length(roots)),vd.dim,vd.dim)
        
        vectors = roots[labels]
        if rank(matrix(vd.field,hcat(vectors...))) == vd.dim
            
            gram = Gram[labels,labels]
            gram_hash = Gram_hash[labels,labels]
            
            simple_graph_plus_colors = to_SimpleGraph_plus_colors(Gram,labels)
               
            group = sort([collect(sort(g)) for g in eachcol(gram_hash)]) |> collect
            if group in keys(grouped_candidates)
                push!(grouped_candidates[group],(labels,gram,simple_graph_plus_colors))
            else
                push!(
                    grouped_candidates,
                    group => [(labels,gram,simple_graph_plus_colors)]
                )
            end
        end


    end

    println("Number of candidate diagrams: $(sum(length.(values(grouped_candidates))))")
    
    id = identity_matrix(vd.field,vd.dim)
    
    for candidates in values(grouped_candidates), pair in Combinatorics.powerset(candidates,1,2)
        
        c1 = pair[1]
        c2 = length(pair) == 2 ? pair[2] : pair[1]

        (labels1,gram1,graph1) = c1
        (labels2,gram2,graph2) = c2
   
        pairings = graph_pairings(vd,graph1,graph2)
        
        for p in pairings
            
            labels_pairs = [(labels1[p],labels2[q]) for (p,q) in p]
                
            m1 = matrix(vd.field,hcat(roots[[p[1] for p in labels_pairs]]...))
            m2 = matrix(vd.field,hcat(roots[[p[2] for p in labels_pairs]]...))


            t = m1 * inv(m2)


            if t≠id && is_integral(vd,t) && preserves_upper_sheet(vd,t) && is_infinite_order_isometry(vd,t,true,true,true,true,true)
                display(p)
                return true
            end


        end
    end

    return false


 end



function inf_ord_sym_n_walls(vd,roots,das,type=:vertices)
   
    Gram = [c' * vd.gram_matrix.entries * d for c in roots, d in roots]
    Gram_hash = hash.(Gram) 


    subsets = begin
        if type === :vertices
            collect.(CoxeterDiagrams.all_spherical_of_rank(das,vd.dim-1))
        elseif type === :n_walls
            Combinatorics.powerset(collect(1:length(roots)),vd.dim-1,vd.dim-1)
        end

    end

    grouped_candidates = Dict{
                              Vector{Vector{UInt64}},
                              Vector{Tuple{Vector{Int64}, Vector{nf_elem}, Matrix{nf_elem}, Tuple{LightGraphs.SimpleGraphs.SimpleGraph{Int64}, Dict{Any, Any}, Dict{Any, Any}}}}
                             }()
    for labels in subsets 
        
        vectors = roots[labels]
        if rank(matrix(vd.field,hcat(vectors...))) == vd.dim-1
            
            gram = Gram[labels,labels]
            gram_hash = Gram_hash[labels,labels]
            last_vector = intersection_vector(vd,vectors)
            type === :vertices && @assert norm_squared(vd,last_vector) < 0
            
            if norm_squared(vd,last_vector) == 0
                continue
            end
            
            simple_graph_plus_colors = to_SimpleGraph_plus_colors(Gram,labels)
               
            group = sort([collect(sort(g)) for g in eachcol(gram_hash)]) |> collect
            if group in keys(grouped_candidates)
                push!(grouped_candidates[group],(labels,last_vector,gram,simple_graph_plus_colors))
            else
                push!(
                    grouped_candidates,
                    group => [(labels,last_vector,gram,simple_graph_plus_colors)]
                )
            end
        end


    end

    println("Number of candidate diagrams: $(sum(length.(values(grouped_candidates))))")
    
    id = identity_matrix(vd.field,vd.dim)
    
    for candidates in values(grouped_candidates), pair in Combinatorics.powerset(candidates,1,2)
        
        c1 = pair[1]
        c2 = length(pair) == 2 ? pair[2] : pair[1]

        (labels1,last_vec1,gram1,graph1) = c1
        (labels2,last_vec2,gram2,graph2) = c2
        
        
        issq,sqrt = issquare(norm_squared(vd,last_vec1)//norm_squared(vd,last_vec2))
        if !issq || sqrt ∉ vd.ring
            continue
        end
        last_vec2 .*= sqrt
        

        if norm_squared(vd,last_vec1) ≠ norm_squared(vd,last_vec2)
            continue
        end

        pairings = graph_pairings(vd,graph1,graph2)
        
        for p in pairings, sg in [+1,-1]
            
            labels_pairs = [(labels1[p],labels2[q]) for (p,q) in p]
                
            m1 = matrix(vd.field,hcat(last_vec1,roots[[p[1] for p in labels_pairs]]...))
            m2 = matrix(vd.field,hcat(sg*last_vec2,roots[[p[2] for p in labels_pairs]]...))

            type === :vertices && @assert norm_squared(vd,sg*last_vec2) <  0
            @assert all(times(vd,sg*last_vec2, r) == 0 for r in roots[[p[2] for p in labels_pairs]])


            t = m1 * inv(m2)


            if t≠id && is_integral(vd,t) && preserves_upper_sheet(vd,t) && is_infinite_order_isometry(vd,t,true,true,true,true,true)
                display(p)
                return true
            end


        end
    end

    return false
end




function intersection_vector(vd,roots)
    
    rk,mat = nullspace(vcat([matrix(vd.field,r' * vd.gram_matrix.entries) for r in roots]...))
    @assert rk == 1
    v = reshape(mat.entries,vd.dim)
    @assert all(times(vd,v,r) == 0 for r in roots)
    if times(vd,v,basepoint(vd)) > 0
        v = -v
    end
    #@assert norm_squared(vd,v) ≤ 0
    return v 
end

function normalize(vd,vec)
    vec2 =  lcm_denominators(vd.ring,vec) .* vec
    vec2 .//= vd.field(elem_gcd(vd.ring,vec2))
    @assert all(v in vd.ring for v in vec2)
    @assert isunit(vd.ring(elem_gcd(vd.ring,vec2)))
    return vec2
end

function infinite_order_symmetry(vd,roots,das)
    # Assume that we're not working in ℚ: so the vertices are never ideal vertices
    vertices_diagrams = CoxeterDiagrams.all_spherical_of_rank(das,das.d)
    vertices_roots = [roots[collect(vd)] for vd in vertices_diagrams]
    
    
    vr1 = vertices_roots[1] # we fix a vertex that we take as base
    [pairings(vd,vr1,vr2) |> display for vr2 in vertices_roots]
end

function to_SimpleGraph_plus_colors(Gram,labels)

    d = length(labels)

    g = LightGraphs.SimpleGraph(d)
    e = LightGraphs.edgetype(g)
    edges_color = Dict()
    vertices_color = Dict()
    for i in 1:d
        push!(vertices_color,i=>Gram[labels[i],labels[i]])
        for j in i+1:d
            LightGraphs.add_edge!(g,i,j)
            # it seems LightGraphs's algorithm wants both direction for each edge…
            push!(edges_color,e(i,j)=>Gram[labels[i],labels[j]])
            push!(edges_color,e(j,i)=>Gram[labels[j],labels[i]])
        end
    end
    return g, edges_color,vertices_color 
end

function graph_pairings(vd,graph1,graph2)

    g1,ce1,cv1 = graph1 
    g2,ce2,cv2 = graph2 

    ec_rel(e1,e2) = ce1[e1] == ce2[e2]
    vc_rel(v1,v2) = cv1[v1] == cv2[v2]

    return LightGraphs.Experimental.all_isomorph(g1,g2, LightGraphs.Experimental.VF2(),edge_relation=ec_rel,vertex_relation=vc_rel)
    
end

function pairing_to_labels(vd,vr1,vr2,pairing)
    return [(vr1[p][1],vr2[q][1]) for (p,q) in pairing]  
end

function pairing_to_matrix(vd,vr1,vr2,pairing)

    m1 = matrix(vd.field,hcat([x[2] for x in vr1[[p[1] for p in pairing]]]...))
    m2 = matrix(vd.field,hcat([x[2] for x in vr2[[p[2] for p in pairing]]]...))

    return m1 * inv(m2)

end

function is_integral(vd,mat)
    all(coeff ∈ vd.ring for coeff in mat)
end

function preserves_the_form(vd,mat)
    mat' * vd.gram_matrix * mat == vd.gram_matrix
end

function is_integrally_invertible(vd,mat)
    
    if !isinvertible(mat)
        return false
    end

    return is_integral(vd,inv(mat))
end

function preserves_upper_sheet(vd,mat)
    to_diag_rep(vd,mat.entries * basepoint(vd))[1] > 0
end
