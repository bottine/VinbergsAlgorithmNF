

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

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::nf_elem) = Base.isless(x.elem_in_nf,y) 
Base.isless(x::nf_elem, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(x,y.elem_in_nf)

Base.isless(x::NfAbsOrdElem{AnticNumberField,nf_elem}, y::NfAbsOrdElem{AnticNumberField,nf_elem}) = Base.isless(x.elem_in_nf,y.elem_in_nf) 

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

lcm_denominators(ring,iter) = elem_lcm(ring,[Hecke.denominator(x,ring) for x in iter]) 

divides(a::NfAbsOrdElem{AnticNumberField,nf_elem},b::NfAbsOrdElem{AnticNumberField,nf_elem},ring) = b.elem_in_nf//a.elem_in_nf ∈ ring
divides(a,b,ring) = b//a∈ring

function t2_exact(x::S) where S <: NumFieldElem
    @assert istotally_real(parent(x))
    return trace(x^2)
end
function t2_exact(x::NfAbsOrdElem)
  return t2_exact(x.elem_in_nf)
end

# thanks Tommy Hofmann
function dot(basis::Vector{NfOrdElem},coefficients::Vector{Int64},tmp::nf_elem,dest=basis[1].parent.nf(0))
    for (b,c) in zip(basis,coefficients)
        mul!(tmp,b.elem_in_nf,c)
        add!(dest,dest,tmp)
    end
    return dest
end

function short_t2_elems(O::NfAbsOrd, lb, ub)
    

    @toggled_assert istotally_real(nf(O))

    trace = Hecke.trace_matrix(O)
    basis = Hecke.basis(O)
    
    tmp = O.nf(0)
    lat = Hecke.short_vectors(Zlattice(gram = trace), lb, ub)
    candidates = [(VinbergsAlgorithmNF.dot(basis,v,tmp),t) for (v,t) in lat]

    @toggled_assert all(t2_exact(c)==t for (c,t) in candidates)
    @toggled_assert all(lb-1 ≤ t2(c) && t2(c) ≤ ub+1 for (c,t) in candidates)
    return candidates
end



function non_neg_short_t2_elems(O::NfAbsOrd, lb, ub)
    candidates = short_t2_elems(O,lb,ub)
    if lb ≤ 0
        push!(candidates,(O.nf(0),fmpq(0)))
    end
    
    return map(v -> (abs(v[1]),v[2]), candidates)
end

approx(x, abs_tol::Int = 32) = BigFloat(conjugates_real(x,abs_tol)[1])
approx(x::NfAbsOrdElem{AnticNumberField,nf_elem}, abs_tol::Int = 32) = BigFloat(conjugates_real(x.elem_in_nf,abs_tol)[1])

approx_sum_at_places(val;first_place_idx) = sum(convert.(Float64,conjugates_real(val,32)[first_place_idx:end]))

# Thanks Tommy Hofmann
function get_enclosing_interval(x::arb)
    a, b = BigFloat(), BigFloat()
    ccall((:arb_get_interval_mpfr, Hecke.libarb), Cvoid, (Ref{BigFloat}, Ref{BigFloat}, Ref{arb}), a, b, x)
    return a, b
end

get_enclosing_intervals_float64(x::nf_elem, abs_tol::Int = 64) = [Float64(lower,RoundDown)..Float64(upper,RoundUp) for (lower, upper) in get_enclosing_interval.(conjugates_real(x,abs_tol))]::Vector{IntervalArithmetic.Interval{Float64}}
lower_float64(x::arb) = Float64(get_enclosing_interval(x)[1],RoundDown)
upper_float64(x::arb) = Float64(get_enclosing_interval(x)[2],RoundUp)

lower_float64(x::nf_elem) = lower_float64(conjugates_real(x,64)[1]) 
upper_float64(x::nf_elem) = upper_float64(conjugates_real(x,64)[1]) 

mutable struct boxed_nf_elem
    x::nf_elem
    box::IntervalArithmetic.Interval{Float64}
end

box(x::nf_elem) = boxed_nf_elem(x,get_enclosing_intervals_float64(x)[1])
approx(x::boxed_nf_elem) = x.box

# exact
ex(x::boxed_nf_elem) = x.x
# approximate
ap(x::boxed_nf_elem) = x.box
lo(x::boxed_nf_elem) = x.box.lo
hi(x::boxed_nf_elem) = x.box.hi

Base.:(+)(x::boxed_nf_elem,y::boxed_nf_elem) = boxed_nf_elem(x.x+y.x, x.box + y.box)
Base.:(*)(x::boxed_nf_elem,y::boxed_nf_elem) = boxed_nf_elem(x.x*y.x, x.box * y.box)
Base.:(-)(x::boxed_nf_elem,y::boxed_nf_elem) = boxed_nf_elem(x.x-y.x, x.box - y.box)
Base.:(-)(x::boxed_nf_elem) = boxed_nf_elem(-x.x, -x.box)
Base.:(//)(x::boxed_nf_elem,y::boxed_nf_elem) = boxed_nf_elem(x.x//y.x, x.box / y.box)

Base.:(==)(x::boxed_nf_elem,y::boxed_nf_elem) = ex(x) == ex(y)
Base.:(==)(x::boxed_nf_elem,y::nf_elem) = ex(x) == y
Base.:(==)(x::boxed_nf_elem,y::Integer) = ex(x) == y

# means that the boxes don't intersect and are on the correct side of each other
Base.:(<<)(x::boxed_nf_elem,y::boxed_nf_elem) = hi(x) < lo(y) 
Base.:(>>)(x::boxed_nf_elem,y::boxed_nf_elem) = hi(y) < lo(x) 

function diagm(K::AnticNumberField,diag)
    n = length(diag)
    M = fill(K(0),n,n)
    for i in 1:n
        M[i,i] = K(diag[i])
    end
    return K.(M)
end

function diagm(diag::Vector{nf_elem})
    @assert length(diag) > 0 "Need a non empty vector"

    return diagm(parent(diag[1]),diag)
end
