function inf_ord_sym2(vd,roots,das)
    
    candidates = Combinatorics.powerset(roots,vd.dim,vd.dim) |> collect
   
    #display(candidates)

    # candidate sets of vectors basis R^{n+1}
    filter!(
        x->isinvertible(matrix(vd.field,hcat(x...))),
        candidates
    )

    for pair in Combinatorics.powerset(candidates,1,2)
        
        c1 = pair[1]
        c2 = length(pair) == 2 ? pair[2] : pair[1]

        # check that c1 and c2 don't share all vectors?

        transfos = pairings(vd,c1,c2) 
        for t in transfos
            
            if is_integral(vd,t) && t≠t^2 # t≠t^2 is a stupid shortcut to check that t is not the identity
                s = inv(t)
                @assert is_integral(vd,s)
                @assert preserve_the_form(vd,t)
                
                # Thanks Tommy Hofmann
                L = splitting_field(minpoly(t))
                tt = change_base_ring(L, t)
                jnf_tt, change = jordan_normal_form(tt)
                if !isdiagonal(jnf_tt)
                    display(jnf_tt) # not diagonal => t is not diagonalizable => infinite order
                end

            else
            end
        end
    end

end

function fixed_points_intersection(vd,tt)
    I_ = matrix(vd.field,LinearAlgebra.I(vd.dim))
    n,N = nullspace(vcat([t-I_ for t in tt]))
    return n,N
end

function infinite_order_symmetry(vd,roots,das)
    # Assume that we're not working in ℚ: so the vertices are never ideal vertices
    vertices_diagrams = CoxeterDiagrams.all_spherical_of_rank(das,das.d)
    vertices_roots = [roots[collect(vd)] for vd in vertices_diagrams]
    
    
    vr1 = vertices_roots[1] # we fix a vertex that we take as base
    [pairings(vd,vr1,vr2) |> display for vr2 in vertices_roots]
end

function vertex_diagram_intersection(vd,roots)
    Hecke._inter([vd.gram*r for r in roots]) 
end


function pairings(vd,vr1,vr2)
    
    function to_SimpleGraph_plus_colors(vertex_roots)

        d = length(vertex_roots)

        g = LightGraphs.SimpleGraph(d)
        e = LightGraphs.edgetype(g)
        edges_color = Dict()
        vertices_color = Dict()
        for i in 1:d
            push!(vertices_color,i=>norm_squared(vd,vertex_roots[i]))
            for j in i+1:d
                LightGraphs.add_edge!(g,i,j)
                # it seems LightGraphs's algorithm wants both direction for each edge…
                push!(edges_color,e(i,j)=>Gram_coeff(vd,vertex_roots[i],vertex_roots[j]))
                push!(edges_color,e(j,i)=>Gram_coeff(vd,vertex_roots[i],vertex_roots[j]))
            end
        end
        return g, edges_color,vertices_color 
    end

    g1,ce1,cv1 = to_SimpleGraph_plus_colors(vr1)
    g2,ce2,cv2 = to_SimpleGraph_plus_colors(vr2)

    ec_rel(e1,e2) = ce1[e1] == ce2[e2]
    vc_rel(v1,v2) = cv1[v1] == cv2[v2]

    pairings =  LightGraphs.Experimental.all_isomorph(g1,g2, LightGraphs.Experimental.VF2(),edge_relation=ec_rel,vertex_relation=vc_rel)
    
    matrices = [pairing_to_matrix(vd,vr1,vr2,p) for p in pairings]
    
    return matrices

end


function pairing_to_matrix(vd,vr1,vr2,pairing)

    m1 = matrix(vd.field,hcat(vr1[[p[1] for p in pairing]]...))
    m2 = matrix(vd.field,hcat(vr2[[p[2] for p in pairing]]...))

    return m1 * inv(m2)

end

function is_integral(vd,mat)
    all(coeff ∈ vd.ring for coeff in mat)
end

function preserve_the_form(vd,mat)
    mat' * vd.gram_matrix * mat == vd.gram_matrix
end

function is_integrally_inversible(vd,mat)
    inv_mat = inv(mat)
    return is_integral(inv_mat)
end
