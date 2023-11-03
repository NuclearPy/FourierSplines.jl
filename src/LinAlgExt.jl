mutable struct SplineLinOp{T<:Number}
    knots::Knots
    coefficients::Array{T}
end

function Base.:+(M1::SplineLinOp, M2::SplineLinOp)
    K = M1.knots
    if K == M2.knots
        return SplineLinOp(K, M1+M2)
    else
        return error("Knots must be the same")
    end
end

function Base.:-(M1::SplineLinOp, M2::SplineLinOp)
    K = M1.knots
    if K == M2.knots
        return SplineLinOp(K, M1-M2)
    else
        return error("Knots must be the same")
    end
end

Base.:*(α::Number, M::SplineLinOp) = SplineLinOp(M.knots, α*M.coefficients)
Base.:*(M::SplineLinOp, α::Number) = SplineLinOp(M.knots, α*M.coefficients)
Base.:/(M::SplineLinOp, α::Number) = SplineLinOp(M.knots, M.coefficients/α)

function Base.:*(M::SplineLinOp, v::Spline)
    K = M.knots
    if K == v.knots
        return SplineLinOp(K, M.coefficients*v.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:*(M1::SplineLinOp, M2::SplineLinOp)
    K = M1.knots
    if K == M2.knots
        return SplineLinOp(K, M1.coefficients*M2.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:\(M::SplineLinOp, v::Spline)
    K = M.knots
    if K == v.knots
        return SplineLinOp(K, M.coefficients\v.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:\(M1::SplineLinOp, M2::SplineLinOp)
    K = M1.knots
    if K == M2.knots
        return SplineLinOp(K, M1.coefficients\M2.coefficients)
    else
        return error("Knots must be the same")
    end
end

for f ∈ [:inv, :exp, :sin, :cos, :tan]
    @eval LinearAlgebra.$f(M::SplineLinOp) = SplineLinOp(M.knots, $f(M.coefficients))
end

LinearAlgebra.det(M::SplineLinOp) = det(M.coefficients)

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