
module Lat
    
    """
        ⊕(M₁,M₂)
    Compute the "block diagonal" matrix with blocks ``M₁`` and ``M₂``, assuming both are square symetric.
    """
    function ⊕(M₁,M₂)

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

    U() = [0 1; 1 0]
    A(n) = Matrix(LinearAlgebra.SymTridiagonal([2 for i in 1:n], [-1 for i in 1:n-1]))
    I(n) = Matrix(LinearAlgebra.I(n))
    B(n) = begin
        @assert n ≥ 3
        M = A_(n)
        @assert false "TODO I'm not sure what the matrix should be!!"
    end
    D(n) = begin
        @assert n≥4
        M = A_(n)
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
# Can't do this here: https://github.com/Nemocas/Nemo.jl/issues/810
#Krat,arat = Hecke.rationals_as_number_field()
#K2,a2 = Hecke.quadratic_field(2)
#K5,a5 = Hecke.quadratic_field(5)

#Bug8 = (a5-1) .* K5.(I_)  ⊕ K5.(gram_E8)


