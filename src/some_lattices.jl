
module Lat
    
    using LinearAlgebra



    """
        ⊕(M₁,M₂)
    Compute the "block diagonal" matrix with blocks ``M₁`` and ``M₂``, assuming both are square symmetric.
    """
    function ⊕(M₁,M₂)

        function good_form(M::Matrix)
            @assert permutedims(M) == M "Only accept symmetric matrices"
            return M
        end

        function good_form(D::Vector)
            M = zeros(Int,length(D),length(D))
            for (i,d) in enumerate(D)
                M[i,i] = d
            end
            return M
        end

        function good_form(c)
            return Matrix(reshape([c],1,1))
        end

        M₁ = good_form(M₁)
        M₂ = good_form(M₂)
        
        n₁ = size(M₁)[1]
        n₂ = size(M₂)[1]

        P = Array{Any, 2}(undef, n₁+n₂, n₁+n₂)
        P[:,:] = fill(0,n₁+n₂,n₁+n₂)
        P[1:n₁,1:n₁] = M₁
        P[n₁+1:end,n₁+1:end] = M₂

        return P
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


