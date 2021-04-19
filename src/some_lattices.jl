
module Lat
    
    using LinearAlgebra

    """
        ⊕(M₁,M₂)
    Compute the "block diagonal" matrix with blocks ``M₁`` and ``M₂``, assuming both are square symetric.
    """
    function ⊕(M₁,M₂)

        if typeof(M₁) <: Vector
            @assert length(M₁) == 1
            M₁ = Matrix(reshape(M₁,1,1))
        end
        if !(typeof(M₁) <: Matrix)
            M₁ = Matrix(reshape([M₁],1,1))
        end
        if typeof(M₂) <: Vector
            @assert length(M₂) == 1
            M₂ = Matrix(reshape(M₂,1,1))
        end
        if !(typeof(M₂) <: Matrix)
            M₂ = Matrix(reshape([M₂],1,1))
        end

        n₁ = size(M₁)[1]
        n₂ = size(M₂)[1]
        @assert size(M₁) == (n₁,n₁)
        @assert size(M₂) == (n₂,n₂)
        @assert M₁' == M₁
        @assert M₂' == M₂

        if  n₂ == 0
            return M₁
        elseif n₁ == 0
            return M₂
        else
            P = Array{Any, 2}(undef, n₁+n₂, n₁+n₂)
            P[:,:] = fill(0,n₁+n₂,n₁+n₂)
            P[1:n₁,1:n₁] = M₁
            P[n₁+1:end,n₁+1:end] = M₂

            return P
        end
    end

    block_diag(M...) = reduce(⊕,M)  

    U() = [0 1; 1 0]
    A(n) = Matrix(LinearAlgebra.SymTridiagonal([2 for i in 1:n], [-1 for i in 1:n-1]))
    I(n) = Matrix(Int.(LinearAlgebra.I(n)))
    B(n) = begin
        @assert n ≥ 3
        M = A(n)
        M[end,end] = 1
        M[end-1,end] = -1
        M[end,end-1] = -1
        return M
        #@assert false "TODO I'm not sure what the matrix should be!!"
    end
    D(n) = begin
        @assert n≥4
        M = A(n)
        M[end,end-1] = 0
        M[end-1,end] = 0
        M[end,end-2] = -1
        M[end-2,end] = -1
        return M
    end

    E(n) = begin
        if n==6
            return [ 
                2  0 -1  0  0  0;
                0  2  0 -1  0  0;
                -1  0  2 -1  0  0;
                0 -1 -1  2 -1  0;
                0  0  0 -1  2 -1;
                0  0  0  0 -1  2
            ]
        elseif n == 7
            return [
                2  0 -1  0  0  0  0;
                 0  2  0 -1  0  0  0;
                -1  0  2 -1  0  0  0;
                 0 -1 -1  2 -1  0  0;
                 0  0  0 -1  2 -1  0;
                 0  0  0  0 -1  2 -1;
                 0  0  0  0  0 -1  2
            ]
        elseif n == 8
            return [ 
            2  0 -1  0  0  0  0  0;
             0  2  0 -1  0  0  0  0;
            -1  0  2 -1  0  0  0  0;
             0 -1 -1  2 -1  0  0  0;
             0  0  0 -1  2 -1  0  0;
             0  0  0  0 -1  2 -1  0;
             0  0  0  0  0 -1  2 -1;
             0  0  0  0  0  0 -1  2
            ] 
        end
        @assert false "Invalid number"

    end


end


