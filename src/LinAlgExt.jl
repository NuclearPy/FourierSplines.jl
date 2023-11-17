mutable struct SplineLinOp{T<:Number}
    domain::Knots
    range::Knots
    coefficients::Array{T}
end

# Utilities and basic operators

# Logic
Base.:(==)(M1::SplineLinOp, M2::SplineLinOp) = (M1.domain == M2.domain) && (M1.coefficients == M2.coefficients)

# Initialization
Base.zeros(T::Type, Domain::Knots, Range::Knots) = SplineLinOp(Domain, Range, zeros(T, dimension(Range), dimension(Domain)))
Base.zeros(T::Type, Domain::Knots, Range::Knots) = SplineLinOp(Domain, Range, ones(T, dimension(Range), dimension(Domain)))
Base.fill(value::Number, Domain::Knots, Range::Knots) = SplineLinOp(Domain, Range, fill(value, dimension(Range), dimension(Domain)))

# Copying
Base.copy(M::SplineLinOp) = SplineLinOp(M.domain, M.range, copy(M.coefficients))

# Indexing, sizes, etc
#Base.parent(S::Spline) = S.coefficients
#Base.parent(K::Knots) = K.knots
#Base.length(S::Spline) = length(S.coefficients)
#Base.size(S::Spline) = size(S.coefficients)
#Base.getindex(S::Spline, i::Int) = S.coefficients[i]
#Base.setindex!(S::Spline, val::Number, i::Int) = (S.coefficients[i] = val)
#Base.axes(S::Spline) = 1:length(S.coefficients)

function Base.:+(M1::SplineLinOp, M2::SplineLinOp)
    K1 = M1.domain
    K2 = M1.range
    if (K1 == M2.domain) && (K2 == M2.range)
        return SplineLinOp(K1, K2, M1+M2)
    else
        return error("Knots must be the same")
    end
end

function Base.:-(M1::SplineLinOp, M2::SplineLinOp)
    K1 = M1.domain
    K2 = M1.range
    if (K1 == M2.domain) && (K2 == M2.range)
        return SplineLinOp(K1, K2, M1-M2)
    else
        return error("Knots must be the same")
    end
end

Base.:*(α::Number, M::SplineLinOp) = SplineLinOp(M.domain, M.range, α*M.coefficients)
Base.:*(M::SplineLinOp, α::Number) = SplineLinOp(M.domain, M.range, α*M.coefficients)
Base.:/(M::SplineLinOp, α::Number) = SplineLinOp(M.domain, M.range, M.coefficients/α)

function Base.:*(M::SplineLinOp, v::Spline)
    K = M.domain
    if K == v.knots
        return Spline(M.range, M.coefficients*v.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:*(M1::SplineLinOp, M2::SplineLinOp)
    K = M1.domain
    if K == M2.range
        return SplineLinOp(M2.domain, M1.range, M1.coefficients*M2.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:\(M::SplineLinOp, v::Spline)
    K = M.range
    if K == v.knots
        return Spline(M.domain, M.coefficients\v.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:\(M1::SplineLinOp, M2::SplineLinOp)
    K = M1.range
    if K == M2.range
        return SplineLinOp(M2.domain, M1.domain, M1.coefficients\M2.coefficients)
    else
        return error("Knots must be the same")
    end
end

Base.adjoint(M::SplineLinOp) = SplineLinOp(M.range, M.domain, parent(M'))

Base.:/(M1::SplineLinOp, M2::SplineLinOp) = (M2' \ M1')'

Base.:^(M::SplineLinOp, α::Number) = SplineLinOp(M.domain, M.range, M.coefficients^α)

for f ∈ [:inv, :exp, :sin, :cos, :tan, :sqrt]
    @eval LinearAlgebra.$f(M::SplineLinOp) = SplineLinOp(M.domain, M.range, $f(M.coefficients))
end

LinearAlgebra.det(M::SplineLinOp) = det(M.coefficients)
LinearAlgebra.logdet(M::SplineLinOp) = logdet(M.coefficients)

# Norms and operator norms

#Spline and SplineLinOp
LinearAlgebra.norm(S::Spline) = LinearAlgebra.norm(S.coefficients, Inf)
LinearAlgebra.norm(S::Spline, p::Real) = LinearAlgebra.norm(S.coefficients, p)

LinearAlgebra.norm(P::FourierSpline) = sum(LinearAlgebra.norm(P[:,n]) for n in axes(P)[2])

function LinearAlgebra.norm(P::FourierSpline, p::Real)
    N = lastindex(P)
    v = zeros(2*N+1)
    for n in 1:(2*N+1)
        v[n] = LinearAlgebra.norm(P[:, n-N-1])
    end
    return LinearAlgebra.norm(v, p)
end

function LinearAlgebra.norm(P::FourierSpline, p::Real, q::Real)
    N = lastindex(P)
    v = zeros(2*N+1)
    for n in 1:(2*N+1)
        v[n] = LinearAlgebra.norm(P[:, n-N-1], q)
    end
    return LinearAlgebra.norm(v, p)
end