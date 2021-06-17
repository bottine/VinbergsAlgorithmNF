function inf_ord_sym2(vd,roots,das)
    
    candidates = filter(x -> length(Set(x))==vd.dim, Base.product([roots for i in 1:vd.dim]...) |> collect)

    # candidate sets of vectors basis R^{n+1}
    filter!(
        x->det(x)≠0,
        candidates
    )

    display(candidates[1])

    for c1 in candidates, c2 in candidates

        transfos = pairings(vd,c1,c2)  
        
        for t in transfos
            
            if is_integral(vd,t)
                # Here means: t is invertible of inverse s with is_integral(vd,s) holding
                # and t preserves the form since t sends c1 to c2 and preserve the Gram matrices

                if any(abs(t) ≠1 for t in eigen(t))
                    display(c1)
                    display(c2)
                    display(t)
                    println("------------------------") 
                end
                
            end


        end

    end

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

    @assert issorted([p[1] for p in pairing])


    m1 = reduce(vcat,vr1)
    display(m1)
    m2 = reduce(vcat,vr2[[p[2] for p in pairing]])

    return m1 * inv(m2)

end

function is_integral(vd,mat)
    all(coeff ∈ vd.ring for coeff in mat)
end

function preserve_the_form(vd,mat)
    mat * vd.gram * mat' == vd.gram
end

function is_integrally_inversible(vd,mat)
    inv_mat = inv(mat)
    return is_integral(inv_mat)
end
