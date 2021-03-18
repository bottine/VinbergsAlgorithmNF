
function diagonalize(ring,A::Matrix)
    # returns T and D with D = T'GT
    # algorithm copied from there https://math.stackexchange.com/questions/1388421/reference-for-linear-algebra-books-that-teach-reverse-hermite-method-for-symmetr
    # plus a gcd step to reduce the growth of values
    # plus the "Case 0" step to further reduce, but it's not enough
    @assert LinearAlgebra.issymmetric(A) "A must be symmetric."



    (n,m) = size(A)
    @assert n == m

    i0 = 1
    I_ = LinearAlgebra.I 
    M = ring.([A I_])
    while i0 ≤ n
       
      
        # look at non zero diagonal entries
        non_zero_diag = [k for k in i0:n if M[k,k] ≠ ring(0)]
        non_zero_diag = sort!(non_zero_diag,by=(k -> abs(M[k,k])))

#        println("====================")
#        println("i0 = $i0")
#        display(M)
#        println("")
#        println("non_zero_diag = $non_zero_diag")
#    
        if length(non_zero_diag) == 0
            non_zero_coordinates = [(i,j) for  i in i0:n, j in i0:n if M[i,j]≠0]
            if isempty(non_zero_coordinates)
                break
            end
            (i,j) = (sort(non_zero_coordinates, by=(x-> abs(M[x[1],x[2]]))))[1]
            M[i,:] = M[i,:] + M[j,:]
            M[:,i] = M[:,i] + M[:,j]
        else
            
            k = non_zero_diag[1]
            M[i0,:], M[k,:] = M[k,:], M[i0,:]
            M[:,i0], M[:,k] = M[:,k], M[:,i0]

            for i in i0+1:n
                g =  elem_gcd(ring,[M[i0,i0],M[i0,i]])
                mizi = ring(M[i0,i].elem_in_nf//g.elem_in_nf)
                miziz = ring(M[i0,i0].elem_in_nf//g.elem_in_nf)
                M[i,:] = (-mizi .* M[i0,:] + miziz .* M[i,:])
                M[:,i] = (-mizi .* M[:,i0] + miziz .* M[:,i])
            end
            i0 = i0 + 1
        end
    end
   

    D = M[1:n,1:n]
    Q = M[1:n,n+1:2*n]
    P = Q'
   

    @assert LinearAlgebra.isdiag(D) "D is diagonal", D
    @assert P'*A*P == D "We have a diagonalization (part 1)"
    #@assert A == inv(P')*D*inv(P) "We have a diagonalization (part 2)"
    
    return (D,P)

end

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

