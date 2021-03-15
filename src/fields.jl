

function Base.:>(x::nf_elem,y::nf_elem,p::InfPlc)
    @assert parent(x) == parent(y) && istotally_real(parent(x)) ""
    
    x == y && return false
    return ispositive(x-y,p)
end
function Base.:≥(x::nf_elem,y::nf_elem,p::InfPlc)
    @assert parent(x) == parent(y) && istotally_real(parent(x)) "Need totally real common parent field"
    
    x == y && return true
    return ispositive(x-y,p)
end

Base.:>(x::nf_elem,y,p::InfPlc) = >(x,parent(x)(y),p)
Base.:>(x,y::nf_elem,p::InfPlc) = >(parent(y)(x),y,p)
Base.:≥(x::nf_elem,y,p::InfPlc) = ≥(x,parent(x)(y),p)
Base.:≥(x,y::nf_elem,p::InfPlc) = ≥(parent(y)(x),y,p)

Base.:<(x,y,p::InfPlc) = >(y,x,p)
Base.:≤(x,y,p::InfPlc) = ≥(y,x,p)

function Base.isless(x::nf_elem,y::nf_elem)
    places = infinite_places(parent(x))
    return <(x,y,places[1])
end

Base.isless(x::nf_elem,y) =  Base.isless(x,parent(x)(y))
Base.isless(x,y::nf_elem) = Base.isless(parent(y)(x),y)

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem},y) =  Base.isless(x,parent(x).nf(y))
Base.isless(x,y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(parent(y).nf(x),y)

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::nf_elem) = Base.isless(parent(x).nf(x),y) 
Base.isless(x::nf_elem, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(x,parent(y).nf(y))

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(parent(x).nf(x),parent(y).nf(y)) 

Base.abs(x) = x<0 ? -x : x


function ideal_gcd(ring,elems)
    idls = [ideal(ring,ring(e)) for e in elems]
    gcd_idls = reduce(gcd,idls)
    return gcd_idls
end

function elem_gcd(ring,elems)
    is_principal,gcd_elems = Hecke.isprincipal(ideal_gcd(ring,elems))
    @assert is_principal 
    return gcd_elems
    
end

function ideal_lcm(ring,elems)
    idls = [ideal(ring,ring(e)) for e in elems]
    lcm_idls = reduce(lcm,idls)
    return lcm_idls
end

function elem_lcm(ring,elems)
    is_principal,lcm_elems = Hecke.isprincipal(ideal_lcm(ring,elems))
    @assert is_principal 
    return lcm_elems
    
end

divides(a::NfAbsOrdElem{AnticNumberField,nf_elem},b::NfAbsOrdElem{AnticNumberField,nf_elem},ring) = b.elem_in_nf//a.elem_in_nf ∈ ring
divides(a,b,ring) = divides(ring(a),ring(b),ring)

function t2_exact(x::S) where S <: NumFieldElem
    @assert istotally_real(parent(x))
    return trace(x^2)
end
function t2_exact(x::NfAbsOrdElem)
  return t2_exact(x.elem_in_nf)
end

function short_t2_elems(O::NfAbsOrd, lb, ub)
    @assert istotally_real(nf(O))

    trace = Hecke.trace_matrix(O)
    basis = Hecke.basis(O)

    lat = Hecke.short_vectors(Zlattice(gram = trace), lb, ub)
    candidates = [O(Hecke.dot(basis,v)) for (v,t) in lat]



    @toggled_assert all(lb-1 ≤ t2(c) && t2(c) ≤ ub+1 for c in candidates)
    return candidates
end

function non_neg_short_t2_elems(O::NfAbsOrd, lb, ub)
    candidates = short_t2_elems(O,lb,ub)
    if lb == 0
        push!(candidates,O(0))
    end
    map!(abs, candidates, candidates)
    return candidates
end

approx(x) = Float64(conjugates_real(x)[1])
approx(x::NfAbsOrdElem{AnticNumberField,nf_elem}) = Float64(conjugates_real(x.elem_in_nf)[1])

approx_sum_at_places(val;first_place_idx) = sum(convert.(Float64,conjugates_real(val,32)[first_place_idx:end]))

function diagm(K::AnticNumberField,diag)
    n = length(diag)
    M = fill(K(0),n,n)
    for i in 1:n
        M[i,i] = K(diag[i])
    end
    return matrix(K,M)
end

function diagm(diag::Vector{nf_elem})
    @assert length(diag) > 0 "Need a non empty vector"

    return diagm(parent(diag[1]),diag)
end
