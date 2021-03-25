Qx,x = Hecke.QQ["x"]
f = 8*x^3 + 4*x^2 - 4*x - 1
K,a = Hecke.NumberField(f, "a")
L, mL = simplify(K)
ϕ₂ = inv(mL)(a) # = cos(6π/7)
ϕ₀ = 2ϕ₂^2-1     # = cos(2π/7)
ϕ₁ = 2ϕ₀^2-1     # = cos(4π/7)

# In Bugaenko, it's [ϕ,ϕ₁,ϕ₂]
# To match Guglielmetti we have aϕ₀+bϕ₁+cϕ₂ =  [a//2,b//2,c//2]


lattices = [
    ( 
        "diag_-2ϕ₀,1,1,1",
        VA.diagm(L,[-2ϕ₀,1,1,1]),
        [[0, -1, 1, 0], [0, 0, -1, 1], [0, 0, 0, -1], [4*ϕ₂^2, 4*ϕ₂^2 - 4*ϕ₂ - 1, 0, 0], [12*ϕ₂^2 - 2*ϕ₂ - 1, 8*ϕ₂^2 - 2*ϕ₂ - 1, 8*ϕ₂^2 - 2*ϕ₂ - 1, 8*ϕ₂^2 - 2*ϕ₂ - 1]]
    ),
    ( # too slow 
        "diag_-2ϕ₀,1,1,1,1",
        VA.diagm(L,[-2ϕ₀,1,1,1,1]),
        [
            [0, -1, 1, 0, 0],
            [0, 0, -1, 1, 0],
            [0, 0, 0, -1, 1],
            [0, 0, 0, 0, -1],
            [-3*(2ϕ₂)^2 - 2*(2ϕ₂) + 7, -(2ϕ₂)^2 - (2ϕ₂) + 3, 0, 0, 0],
            [-(2ϕ₂) + 1, 3*(2ϕ₂)^2 + (2ϕ₂) - 6, 3*(2ϕ₂)^2 + (2ϕ₂) - 6, 3*(2ϕ₂)^2 + (2ϕ₂) - 6, 0],
            [4*(2ϕ₂)^2 - 6, -(2ϕ₂)^2 - 2*(2ϕ₂) + 4, -(2ϕ₂)^2 - 2*(2ϕ₂) + 4, 2*(2ϕ₂)^2 - 3, 2*(2ϕ₂)^2 - 3],
            [2*(2ϕ₂)^2 - 2*(2ϕ₂) - 1, 2*(2ϕ₂)^2 - (2ϕ₂) - 2, 2*(2ϕ₂)^2 - (2ϕ₂) - 2, 2*(2ϕ₂)^2 - 3, 2*(2ϕ₂)^2 - 3]
        ]
    ),
    (
        "diag_-2ϕ₀,1,etc", # Found in Guglielmetti page 132 WARNING: MY RESULT DOES NOT MATCH HIS
        VA.diagm(L,[-2ϕ₀,1,-4ϕ₀-10ϕ₁-14ϕ₂]),
        [[0, -1, 0],[0, 0, 2*(2ϕ₂)^2 + (2*ϕ₂) - 5],[2*(2ϕ₂)^2 + 2ϕ₂ - 4, (2ϕ₂)^2 - 2, 0],[-2*(2ϕ₂)^2 - 4*ϕ₂ + 5, 0, -3*(2ϕ₂)^2 - 4*ϕ₂ + 7]]
        # For reference, AlVin yields, if I didn't messu up the translation
        # [[0,-1,0],[0,0,-1],[1,-2ϕ₂,0],[-2ϕ₁-4ϕ₂,-2ϕ₂,2ϕ₀],[-2ϕ₀-4ϕ₁-8ϕ₂,0,-2ϕ₁-2ϕ₂]]
    )
]


@testset "Some reflective lattices over ℚ(cos(2π/7))" begin
    
    for (name,matrix,known_roots) in lattices
        @testset "$name" begin
            compare_vinberg_outputs(L,matrix,known_roots)
        end
    end

end

