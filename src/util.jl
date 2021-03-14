
"""
    product(d)

Given the dictionary `d` with entries of the form `a => n_a` with `n_a` a positive integer, compute and return all products of the form ``∏_{a} a^{d_a}`` with ``a`` iterating over the keys of `d` and `0≤d_a≤n_a`.
"""
function products(d)
    ks = collect(keys(d))
    vs = collect(values(d))
    partial_products = Base.product([collect(0:v) for v in vs]...)
    return [prod([k^p for (k,p) in zip(ks,pp)]) for pp in partial_products]
end



# TODO: make exact
function is_necessary_halfspace(cone_roots,root) 
   


    float_cone_roots = Vector{Vector{Float64}}([approx.(cone_root) for cone_root in cone_roots])    
    float_root = Vector{Float64}(approx.(root))
    
    n = length(root) 

    # x' * (A * r) ≤ 0 ∀ r
    # (A * r)' * x ≤ 0 ∀ r

    #x = Variable(n, IntVar)
    x = Variable(n)
    p = satisfy()       # satisfiability question 
    for cone_root in float_cone_roots
        p.constraints += x' * cone_root ≤ 0 # hyperplanes defining the cone
    end
    p.constraints += x' * float_root ≥ 1 # other side of the half space defined by root
    # it should only be strictly bigger than zero, but Convex.jl does not do "strictly", so we change it to ≥ 1 (and since we have a cone, it should be the same result)

    
    Convex.solve!(p,Cbc.Optimizer(verbose=0,loglevel=0), verbose=false, warmstart=false)
    #solve!(p,COSMO.Optimizer(verbose=false), verbose=false)
   

    if p.status == MathOptInterface.INFEASIBLE 
        return false
    elseif p.status == MathOptInterface.OPTIMAL
        #println(p.optval)
        return true
    else
        println("can't tell! ($(p.status))")
        println("             $(p))")
    end

end


function drop_redundant_halfspaces(
    roots
) 
    
    for i in length(roots):-1:1
       
        rr = copy(roots)
        r = popat!(rr,i)

        if ! is_necessary_halfspace(rr,r)
            return drop_redundant_halfspaces(rr) 
        end
    end
    
    return roots
    
end

