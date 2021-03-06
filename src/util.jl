

function diagonalize(ring,A::Matrix,v₀=nothing)
    # returns T and D with D = T'GT
    # algorithm copied from there https://math.stackexchange.com/questions/1388421/reference-for-linear-algebra-books-that-teach-reverse-hermite-method-for-symmetr
    # plus a gcd step to reduce the growth of values
    # plus the "Case 0" step to further reduce, but it's not enough
    @assert LinearAlgebra.issymmetric(A) "A must be symmetric."



    (n,m) = size(A)
    @assert n == m

    i0 = 1
    I_ = ring.(collect(LinearAlgebra.I(n)))
    M = ring.([A I_])

    if v₀ ≠ nothing
        @assert v₀' * A * v₀ ≠ 0
        
        # Do dumb stuff to get v₀ as first vector: this can probably be made cleaner
        # We essentially just force v₀ as first vector
        
        k = [i for i in 1:n if v₀[i]≠0][1]::Int

        M[k,:] = v₀[k].*M[k,:]
        M[:,k] = v₀[k].*M[:,k]

        for i in 1:n
            if i ≠ k
                M[k,:] = M[k,:] +  v₀[i] .* M[i,:]
                M[:,k] = M[:,k] +  v₀[i] .* M[:,i]
            end
        end
      
        M[1,:], M[k,:] = M[k,:], M[1,:]
        M[:,1], M[:,k] = M[:,k], M[:,1]
        
        @assert M[1,n+1:end] == v₀
    end

    D = M[1:n,1:n]
    Q = M[1:n,n+1:2*n]
    P = Q'


    @assert P'*A*P == D 
    
    while i0 ≤ n
       
      
        # look at non zero diagonal entries
        non_zero_diag = [k for k in i0:n if M[k,k] ≠ ring(0)]

    
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

function diagonalize_in_field(field,A::Matrix,v₀=nothing)
    # returns T and D with D = T'GT
    # algorithm copied from there https://math.stackexchange.com/questions/1388421/reference-for-linear-algebra-books-that-teach-reverse-hermite-method-for-symmetr
    # plus the "Case 0" step to further reduce, but it's not enough
    @assert LinearAlgebra.issymmetric(A) "A must be symmetric."



    (n,m) = size(A)
    @assert n == m

    i0 = 1
    I_ = field.(collect(LinearAlgebra.I(n)))
    M = field.([A I_])

    if v₀ ≠ nothing
        @assert v₀' * A * v₀ ≠ 0
        
        # Do dumb stuff to get v₀ as first vector: this can probably be made cleaner
        # We essentially just force v₀ as first vector
        
        k = [i for i in 1:n if v₀[i]≠0][1]::Int

        M[k,:] = v₀[k].*M[k,:]
        M[:,k] = v₀[k].*M[:,k]

        for i in 1:n
            if i ≠ k
                M[k,:] = M[k,:] +  v₀[i] .* M[i,:]
                M[:,k] = M[:,k] +  v₀[i] .* M[:,i]
            end
        end
      
        M[1,:], M[k,:] = M[k,:], M[1,:]
        M[:,1], M[:,k] = M[:,k], M[:,1]
        
        @assert M[1,n+1:end] == v₀
    end

    D = M[1:n,1:n]
    Q = M[1:n,n+1:2*n]
    P = Q'


    @assert P'*A*P == D 
    
    while i0 ≤ n
       
      
        # look at non zero diagonal entries
        non_zero_diag = [k for k in i0:n if M[k,k] ≠ 0]

    
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
                M[i,:] = (-M[i0,i] .* M[i0,:] + M[i0,i0] .* M[i,:])
                M[:,i] = (-M[i0,i] .* M[:,i0] + M[i0,i0] .* M[:,i])
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




function diagonalize_and_get_scaling(gram,ring,field,v₀=nothing)

    @assert LinearAlgebra.issymmetric(gram)
    n = size(gram)[1]

    diagonal_values,diagonal_basis = diagonalize(ring,gram,v₀)
    @assert LinearAlgebra.isdiag(diagonal_values)
   
    diagonal_basis_vecs = [[diagonal_basis[i,j] for i in 1:n] for j in 1:n]

    inverse = Hecke.inv(matrix(field,field.(diagonal_basis)))
    scaling = [abs(lcm_denominators(ring,[inverse[i,j] for j in 1:n])) for i in 1:n]
    # clever scaling does not seem to work for now, TODO
    #scaling = [abs(lcm_denominators(ring,inverse)) for i in 1:n]
    return diagonal_basis_vecs, LinearAlgebra.diag(diagonal_values), scaling 
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

⇒(x::Bool,y::Bool) = !x || y
⇔(x::Bool,y::Bool) = (x && y) || (!x && !y)

Base.copy(x::NfAbsOrdElem) = x # Otherwise Nemo throws errors on small matrices



# TODO: make exact
function is_necessary_halfspace(gram,cone_roots,root) 
   
    lp_precision = Options.lp_precision()

    float_cone_roots_grammed = Vector{Vector{BigFloat}}([[approx(r,lp_precision) for r in gram*cone_root] for cone_root in cone_roots])    
    float_root_grammed = Vector{BigFloat}([approx(r,lp_precision) for r in gram*root])

    n = length(root) 

    # x' * (A * r) ≤ 0 ∀ r
    # (A * r)' * x ≤ 0 ∀ r

    #x = Variable(n, IntVar)
    x = Variable(n)
    p = maximize(big(0), numeric_type=BigFloat)       # satisfiability question 
    for float_cone_root_grammed in float_cone_roots_grammed
        p.constraints += x' * float_cone_root_grammed ≤ big(0) # + big(2)^(-256) # hyperplanes defining the cone
    end
    p.constraints += x' * float_root_grammed ≥ big(1) # other side of the half space defined by root
    # it should only be strictly bigger than zero, but Convex.jl does not do "strictly", so we change it to ≥ 1 (and since we have a cone, it should be the same result)

    
    #Convex.solve!(p,Cbc.Optimizer(verbose=0,loglevel=0), verbose=false, warmstart=false)
    #Convex.solve!(p,COSMO.Optimizer(verbose=false), verbose=false, warmstart=false)
    Convex.solve!(p,Tulip.Optimizer{BigFloat}(), verbose=false)
   

    if p.status == MathOptInterface.INFEASIBLE
        return false
    elseif p.status == MathOptInterface.OPTIMAL
        return true
    else
        @error "LP program couldn't decide feasibility."
        println("can't tell! ($(p.status))")
        println("             $(p))")
    end

end


function drop_redundant_halfspaces(
    gram,
    roots
) 
    
    for i in 1:length(roots)
       
        rr = copy(roots)
        r = popat!(rr,i)

        if ! is_necessary_halfspace(gram,rr,r)
            return drop_redundant_halfspaces(gram,rr) 
        end
    end
    
    return roots
    
end



function matrix_to_dot(Mat)
    
    (m,n) = size(Mat)
    
    @assert m == n

    dot_string = """
strict graph {
    layout=neato
    node [shape=point];
    """

    function label_to_edge_type(k)
        if k == 1
            return "[style=dotted,label=$k,weight=0]"
        elseif k == 0
            return "[penwidth=3,label=$k,weight=2]"
        elseif k == 3
            return "[label=$k]"
        elseif k > 3
            return "[color = \"" * "black:invis:"^(k-3) * ":black\",label=$k]"
        end
    end

    for i in 1:m
        for j in i+1:n
            if Mat[i,j] ≠ 2
                dot_string = dot_string * "\t$i -- $j" * label_to_edge_type(Mat[i,j]) * ";\n"
            end
        end
    end

    dot_string = dot_string * "}"


    return dot_string

end

function matrix_to_dot_file(Mat, path)
    s = open(path,"w") do file
        print(file, matrix_to_dot(Mat))
    end
end


##########################################
#
#  Read matrix/lattice info from JSON file
#
##########################################

function JSON_field(desc)
    if desc == "QQ"
        
        K,a = Hecke.rationals_as_number_field()
        return ((K,a),[K(1)])
    
    elseif desc == "RC7"
        Qx,x = Hecke.QQ["x"]
        f = 8*x^3 + 4*x^2 - 4*x - 1
        K,a = Hecke.NumberField(f, "a")
        L, mL = simplify(K)
        ϕ₂ = inv(mL)(a) # = cos(6π/7)
        ϕ₀ = ϕ₂^2-1     # = cos(2π/7)
        ϕ₁ = ϕ₀^2-1     # = cos(4π/7)

        # To match Guglielmetti we have aϕ₀+bϕ₁+cϕ₂ =  [a//2,b//2,c//2]
        
        return ((L,ϕ₂),[1,2ϕ₂,(2ϕ₂)^2-2])

    elseif "Quad" ∈ keys(desc)
        p = desc["Quad"]::Int
        K,a = Hecke.quadratic_field(p)
        b = -a
        ϕ = p % 4 == 1 ? (b+1)//2 : b

        return ((K,b),[K(1),b])
    end
end

function replace_symbol(s, field_gen, ring_gens)
    if s == Symbol("A")
        return field_gen
    end
    for (idx,letter) in enumerate("abcdefghi") # don't need that many letters
        if s == Symbol(letter)
            @assert idx ≤ length(ring_gens) "Can't translate entry: too many letters!!!"
            return ring_gens[idx]
        end
    end
    return s
end

function replace_in_expr!(expr, field_gen, ring_gens)
    expr.head = replace_symbol(expr.head, field_gen, ring_gens)
    for (i,v) in enumerate(expr.args)
        if typeof(v) == Expr
            replace_in_expr!(v,field_gen,ring_gens)
        else
            expr.args[i] = replace_symbol(v, field_gen, ring_gens)
        end
    end
end

function read_matrix_entry(desc, field, field_gen, ring_gens)
    if typeof(desc) <: Integer
        return field(desc)
    elseif typeof(desc) <: AbstractVector
        return sum(ring_gen * field(coeff) for (coeff,ring_gen) in zip(desc,ring_gens))
    elseif typeof(desc) <: AbstractString
        parsed = Meta.parse(desc)
        replace_in_expr!(parsed, field_gen, ring_gens)
        return eval(parsed)
    end
    @assert false "Couldn't parse entry" 
end

function read_lattice_json(desc::String)
    
    j = JSON.parse(desc)

    (field,field_gen),ring_gens = JSON_field(j["field"])
    matrix = vcat([Vector{nf_elem}(map(x -> field(read_matrix_entry(x,field,field_gen,ring_gens)),line))' for line in j["matrix"]]...) 
    reflective = j["reflective"]
    
    if reflective
        Coxeter_matrix = vcat([Vector{Int}(line)' for line in j["matrix"]]...)
        return field,(field_gen,ring_gens),matrix,reflective, Coxeter_matrix
    else
        return field,(field_gen,ring_gens),matrix,reflective
    end
end

function read_lattices_json_file(desc::AbstractString)
    
    jj = JSON.parsefile(desc)
    lattices = []
    for j in jj
        (field,field_gen),ring_gens = JSON_field(j["field"])
        matrix = vcat([Vector{nf_elem}(map(x -> field(read_matrix_entry(x,field,field_gen,ring_gens)),line))' for line in j["matrix"]]...) 
        reflective = j["reflective"]
        
        if reflective
            Coxeter_matrix = vcat([Vector{Int}(line)' for line in j["Coxeter matrix"]]...)
            push!(lattices, (field,(field_gen,ring_gens),matrix,reflective, Coxeter_matrix))
        else
            push!(lattices, (field,(field_gen,ring_gens),matrix,reflective))
        end
    end
    return lattices
end


