module FourierSplines
export Knots, Spline, FourierSpline, knots, order, degree, dimension, coefficients, locateIndex, eval

import LinearAlgebra

# Spline Space
struct Knots{T<:Real}
    knots :: Vector{T}
    order::Int
end

function knots(K::Knots)
    return K.knots
end

function order(K::Knots)
    return K.order
end

function degree(K::Knots)
    return K.order + 1
end

function dimension(K::Knots)
    return Base.length(K.knots)
end

# Logic
Base.:(==)(K1::Knots, K2::Knots) = (K1.knots == K2.knots) && (K1.order == K2.order)

# Spline construction
mutable struct Spline{S <: Real, T <: Number}
    const knots :: Knots{S}
    coefficients :: Vector{T}
end

function knots(S::Spline)
    return knots(S.knots)
end

function coefficients(S::Spline)
    return S.coefficients
end

function order(S::Spline)
    return order(S.knots)
end

function degree(S::Spline)
    return degree(S.knots)
end

LinearAlgebra.norm(S::Spline) = norm(S.coefficients, Inf)
LinearAlgebra.norm(S::Spline, p::Real) = norm(S.coefficients, p)

# Utilities and basic operators

# Logic
Base.:(==)(S1::Spline, S2::Spline) = (S1.knots == S2.knots) && (S1.coefficients == S2.coefficients)

# Initialization
Base.zeros(T::Type, K::Knots) = Spline(K, zeros(T, dimension(K)))
Base.ones(T::Type, K::Knots) = Spline(K, ones(T, dimension(K)))
Base.fill(value::Number, K::Knots) = Spline(K, fill(value, dimension(K)))

# Copying
Base.copy(S::Spline) = Spline(S.knots, copy(S.coefficients))

# Indexing, sizes, etc
Base.parent(S::Spline) = S.coefficients
Base.parent(K::Knots) = K.knots
Base.length(S::Spline) = length(S.coefficients)
Base.size(S::Spline) = size(S.coefficients)
Base.getindex(S::Spline, i::Int) = S.coefficients[i]
Base.setindex!(S::Spline, val::Number, i::Int) = (S.coefficients[i] = val)
Base.axes(S::Spline) = 1:length(S.coefficients)

function locateIndex(S::Spline,t::Real)
    k = knots(S)
    if t<first(k)
        return 0
    elseif t>last(k)
        return length(S)
    else
        for i in 1:(length(S)-1)
            if k[i]<=t && t<=k[i+1]
                return i
            end
        end
    end
end

# Basic operations
function Base.:+(S1::Spline, S2::Spline)
    if S1.knots == S2.knots
        return Spline(S1.knots, S1.coefficients + S2.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:+(S::Spline, α::Number)
    return Spline(S.knots, S.coefficients .+ α)
end

function Base.:+(α::Number, S::Spline)
    return Spline(S.knots, α .+ S.coefficients)
end

function Base.:-(S1::Spline, S2::Spline)
    if S1.knots == S2.knots
        return Spline(S1.knots, S1.coefficients - S2.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:-(S::Spline, α::Number)
    return Spline(S.knots, S.coefficients .- α)
end

function Base.:-(α::Number, S::Spline)
    return Spline(S.knots, α .- S.coefficients)
end

function Base.:*(S1::Spline, S2::Spline)
    if S1.knots == S2.knots
        return Spline(S1.knots, S1.coefficients .* S2.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:*(S::Spline, α::Number)
    return Spline(S.knots, S.coefficients * α)
end

function Base.:*(α::Number, S::Spline)
    return Spline(S.knots, α * S.coefficients)
end

function Base.:/(S1::Spline, S2::Spline)
    if S1.knots == S2.knots
        return Spline(S1.knots, S1.coefficients ./ S2.coefficients)
    else
        return error("Knots must be the same")
    end
end

function Base.:/(S::Spline, α::Number)
    return Spline(S.knots, S.coefficients / α)
end

function Base.:/(α::Number, S::Spline)
    return Spline(S.knots, α ./ S.coefficients)
end

function Base.:^(S::Spline, α::Number)
    return Spline(S.knots, S.coefficients .^ α)
end

function Base.:^(α::Number, S::Spline)
    return Spline(S.knots, α .^ S.coefficients)
end

function Base.:^(S1::Spline, S2::Spline)
    if S1.knots == S2.knots
        return Spline(S1.knots, S1.coefficients .^ S2.coefficients)
    else
        return error("Knots must be the same")
    end
end

# Elementary function evaluations
for f ∈ [:sin, :cos, :tan, :sec, :csc, :cot, :exp, :log, :conj, :sqrt, :abs]
    @eval Base.$f(S::Spline) = Spline(S.knots, $f.(S.coefficients))
end

Base.log(b::Real, S::Spline) = Spline(S.knots, log.(b,S.coefficients))

# Evaluating Splines
function eval(S::Spline, t::Real)
    if order(S) != 1
        return error("Evaluation not yet supported for orders other than 1")
    else
        i=locateIndex(S,t)
        if i == 0 || i == length(S)
            return 0
        else
            K = knots(S)
            t1, t2 = K[i], K[i+1]
            y1, y2=S[i], S[i+1]
            m=(y2-y1)/(t2-t1)
            return m*(t-t1)+y1
        end
    end
end

function (S::Spline)(t::Real) # Calling Spline as a function
    if order(S) != 1
        return error("Evaluation not yet supported for orders other than 1")
    else
        i=locateIndex(S,t)
        if i == 0 || i == length(S)
            return 0
        else
            K = knots(S)
            t1, t2 = K[i], K[i+1]
            y1, y2=S[i], S[i+1]
            m=(y2-y1)/(t2-t1)
            return m*(t-t1)+y1
        end
    end
end

# Fourier-Spline Series
mutable struct FourierSpline{S<:Real, T<:Number}
    knots::Knots{S}
    series::Array{T}
end

# Utilities and basic operators

# Logic
Base.:(==)(P1::FourierSpline, P2::FourierSpline) = (P1.knots == P2.knots) && (P1.series == P2.series)

# Initialization
Base.zeros(T::Type, K::Knots, fourierN::Int) = FourierSpline(K, zeros(T, dimension(K), 2*fourierN+1))
Base.ones(T::Type, K::Knots, fourierN::Int) = FourierSpline(K, ones(T, dimension(K), 2*fourierN+1))
Base.fill(value::Number, K::Knots, fourierN::Int) = FourierSpline(K, fill(value, dimension(K), 2*fourierN+1))

# Copying
Base.copy(P::FourierSpline) = FourierSpline(P.knots, copy(P.series))

# Indexing, sizes, etc
Base.length(P::FourierSpline) = length(P.series)
Base.size(P::FourierSpline) = size(P.series)
Base.lastindex(P::FourierSpline) = trunc(Int, (size(P)[2]-1)/2)
Base.firstindex(P::FourierSpline) = -lastindex(P::FourierSpline)
Base.parent(P::FourierSpline) = P.series
Base.getindex(P::FourierSpline, ::Colon) = P.series[:]

function Base.axes(P::FourierSpline)
    N = lastindex(P)
    M = size(P)[1]
    return (1:M, -N:N)
end

function Base.getindex(P::FourierSpline, i::Int, j::Int)
    N = lastindex(P)
    if -N ≤ j ≤ N
        return P.series[i, j + N + 1]
    else
        a = size(P)
        b = axes(P)
        ind = [i,j]
        return error("BoundsError: Attempt to access ", a, " ", P, " with indices ", b, " at index ", ind)
    end
end

function Base.getindex(P::FourierSpline, ::Colon, j::Int)
    N = lastindex(P)
    if -N ≤ j ≤ N
        return Spline(P.knots, P.series[:, j + N + 1])
    else
        a = size(P)
        b = axes(P)
        ind = [1:dimension(P.knots),j]
        return error("BoundsError: Attempt to access ", a, " ", P, " with indices ", b, " at index ", ind)
    end
end
Base.setindex!(P::FourierSpline, val::Number, i::Int, j::Int) = (P.series[i,j+lastindex(P)+1] = val)
Base.setindex!(P::FourierSpline, val::Spline, ::Colon, j::Int) = (P.series[:,j+lastindex(P)+1] = parent(val))

Base.reshape(P::FourierSpline, sizes...) = reshape(P.series, sizes)

# Basic operations
function Base.:+(P1::FourierSpline, P2::FourierSpline)
    if P1.knots == P2.knots
        N1 = lastindex(P1)
        if N1 ≤ lastindex(P2)
            Sum = copy(P2)
            for n in -N1:N1
                Sum[:,n] += P1[:,n]
            end
            return Sum
        else
            return P2+P1
        end
    else
        return error("Knots must be the same")
    end
end

Base.:*(α::Number, P::FourierSpline) = FourierSpline(P.knots, α*P.series)
Base.:*(P::FourierSpline, α::Number) = FourierSpline(P.knots, α*P.series)

function Base.:-(P1::FourierSpline, P2::FourierSpline)
    if P1.knots == P2.knots
        N2 = lastindex(P2)
        if N2 ≤ lastindex(P1)
            Sum = copy(P1)
            for n in -N2:N2
                Sum[:,n] -= P2[:,n]
            end
            return Sum
        else
            return -1*(P2-P1)
        end
    else
        return error("Knots must be the same")
    end
end

function Base.:*(P1::FourierSpline{T}, P2::FourierSpline{S}) where {T,S}
    K = P1.knots
    R = promote_type(T,S)
    N1 = lastindex(P1)
    N2 = lastindex(P2)
    N = N1 + N2
    if K == P2.knots
        P1P2 = zeros(R, K, N)
        for n in -N:N
            for k in -N1:N1
                if -N2 ≤ n-k ≤ N2
                    P1P2[:,n] += P1[:,k]*P2[:,n-k]
                end
            end
        end
        return P1P2
    else
        return error("Knots must be the same")
    end
end

function truncate(P::FourierSpline{T}, fourierN::Int) where {T}
    if fourierN ≥ lastindex(P)
        return copy(P)
    else
        truncP = zeros(T, P.knots, fourierN)
        for n in -fourierN:fourierN
            truncP[:, n] = copy(P[:,n])
        end
        return truncP
    end
end

# Evaluation
function eval(P::FourierSpline, θ::Real)
    N = lastindex(P)
    fourierExp = 2 * π * θ * im
    val = real.(sum(P[:, n] * exp(fourierExp * n) for n = -N:N))
    return Spline(P.knots, val)
end

function (P::FourierSpline)(θ::Real)
    N = lastindex(P)
    fourierExp = 2 * π * θ * im
    val = real.(sum(P[:, n] * exp(fourierExp * n) for n = -N:N))
    return Spline(P.knots, val)
end

function eval(P::FourierSpline, θ::Real, t::Real)
    return P(θ)(t)
end

function (P::FourierSpline)(θ::Real, t::Real)
    return P(θ)(t)
end

end
