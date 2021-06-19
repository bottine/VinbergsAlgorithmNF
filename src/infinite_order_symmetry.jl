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
)

    # Normally: if it preserve the form and is invertible and integral, the inverse also preserves the form
    # Also, having unit determinant is enough to have an integral inverse
    # But we check too much: better safe than sorry

    check_invertible && @assert isinvertible(t)
    check_integral && @assert is_integral(vd,t)
    check_integrally_invertible && @assert isunit(vd.ring(det(t)))
    check_integrally_invertible && @assert is_integrally_invertible(vd,t)
    check_preserves_form && @assert preserves_the_form(vd,t) 

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

function inf_ord_sym2(vd,roots,das)
    
    candidates = Combinatorics.powerset(collect(enumerate(roots)),vd.dim,vd.dim) |> collect
  
    println("Number of candidate diagrams: $(length(candidates))")
    
    # candidate sets of vectors basis R^{n+1}
    filter!(
        x->isinvertible(matrix(vd.field,hcat([y[2] for y in x]...))),
        candidates
    )

    id = identity_matrix(vd.field,vd.dim)

    println("After dropping non spanning : $(length(candidates))")
    
    for pair in Combinatorics.powerset(candidates,1,2)
        
        c1 = pair[1]
        c2 = length(pair) == 2 ? pair[2] : pair[1]

        # check that c1 and c2 don't share all vectors?

        transfos = pairings(vd,c1,c2) 
        
        for (labels,t) in transfos
           
            if t≠id && is_integral(vd,t) && is_infinite_order_isometry(vd,t,true,true,true,true)
                display(labels)
                return true
            end


        end
    end

    return false
end

function infinite_order_symmetry(vd,roots,das)
    # Assume that we're not working in ℚ: so the vertices are never ideal vertices
    vertices_diagrams = CoxeterDiagrams.all_spherical_of_rank(das,das.d)
    vertices_roots = [roots[collect(vd)] for vd in vertices_diagrams]
    
    
    vr1 = vertices_roots[1] # we fix a vertex that we take as base
    [pairings(vd,vr1,vr2) |> display for vr2 in vertices_roots]
end

function pairings(vd,vr1,vr2)
    
    function to_SimpleGraph_plus_colors(vertex_roots)

        d = length(vertex_roots)

        g = LightGraphs.SimpleGraph(d)
        e = LightGraphs.edgetype(g)
        edges_color = Dict()
        vertices_color = Dict()
        for i in 1:d
            push!(vertices_color,i=>norm_squared(vd,vertex_roots[i][2]))
            for j in i+1:d
                LightGraphs.add_edge!(g,i,j)
                # it seems LightGraphs's algorithm wants both direction for each edge…
                push!(edges_color,e(i,j)=>Gram_coeff(vd,vertex_roots[i][2],vertex_roots[j][2]))
                push!(edges_color,e(j,i)=>Gram_coeff(vd,vertex_roots[i][2],vertex_roots[j][2]))
            end
        end
        return g, edges_color,vertices_color 
    end

    g1,ce1,cv1 = to_SimpleGraph_plus_colors(vr1)
    g2,ce2,cv2 = to_SimpleGraph_plus_colors(vr2)

    ec_rel(e1,e2) = ce1[e1] == ce2[e2]
    vc_rel(v1,v2) = cv1[v1] == cv2[v2]

    pairings =  LightGraphs.Experimental.all_isomorph(g1,g2, LightGraphs.Experimental.VF2(),edge_relation=ec_rel,vertex_relation=vc_rel)
    
    return [(pairing_to_labels(vd,vr1,vr2,p),pairing_to_matrix(vd,vr1,vr2,p)) for p in pairings]

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
