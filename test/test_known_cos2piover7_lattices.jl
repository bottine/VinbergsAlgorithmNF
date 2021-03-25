Qx,x = Hecke.QQ["x"]
f = 8*x^3 + 4*x^2 - 4*x - 1
K,a = Hecke.NumberField(f, "a")
L, mL = simplify(K)
ϕ = inv(mL)(a) # = cos(6π/7)
ψ = 2ϕ^2-1     # = cos(2π/7)
θ = 2ψ^2-1     # = cos(4π/7)

# To match Guglielmetti, the order is thus [ψ,θ,ϕ]


lattices = [
    ( 
        "diag_-2ψ,1,1,1",
        VA.diagm(L,[-2ψ,1,1,1]),
        [[0, -1, 1, 0], [0, 0, -1, 1], [0, 0, 0, -1], [4*ϕ^2, 4*ϕ^2 - 4*ϕ - 1, 0, 0], [12*ϕ^2 - 2*ϕ - 1, 8*ϕ^2 - 2*ϕ - 1, 8*ϕ^2 - 2*ϕ - 1, 8*ϕ^2 - 2*ϕ - 1]]
    ),
    (
        "diag_-2ψ,1,etc", # Found in Guglielmetti page 132 WARNING: MY RESULT DOES NOT MATCH HIS
        VA.diagm(L,[-2ψ,1,-4ψ-10θ-14ϕ]),
        [[0, -1, 0],[0, 0, 2*(2ϕ)^2 + (2*ϕ) - 5],[2*(2ϕ)^2 + 2ϕ - 4, (2ϕ)^2 - 2, 0],[-2*(2ϕ)^2 - 4*ϕ + 5, 0, -3*(2ϕ)^2 - 4*ϕ + 7]]
        # For reference, AlVin yields, if I didn't messu up the translation
        # [[0,-1,0],[0,0,-1],[1,-2ϕ,0],[-2θ-4ϕ,-2ϕ,2ψ],[-2ψ-4θ-8ϕ,0,-2θ-2ϕ]]
    )
]


@testset "Some reflective lattices over ℚ(cos(2π/7))" begin
    
    for (name,matrix,known_roots) in lattices
        @testset "$name" begin
            compare_vinberg_outputs(L,matrix,known_roots)
        end
    end

end

