# Fourier-Spline Series
mutable struct FourierSpline{T<:Number}
    knots::Knots
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

function Base.:+(P::FourierSpline, S::Spline)
    if P.knots == S.knots
        PS = copy(P)
        PS[:,0] += S
        return PS
    else
        return error("Knots must be the same")
    end
end

function Base.:+(S::Spline, P::FourierSpline)
    if P.knots == S.knots
        PS = copy(P)
        PS[:,0] += S
        return PS
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
function evaluate(P::FourierSpline, θ::Real)
    N = lastindex(P)
    fourierExp = 2 * π * θ * im
    val = sum(P[:, n] * exp(fourierExp * n) for n = -N:N)
    return val
end

function (P::FourierSpline)(θ::Real)
    N = lastindex(P)
    fourierExp = 2 * π * θ * im
    val = sum(P[:, n] * exp(fourierExp * n) for n = -N:N)
    return val
end

function evaluate(P::FourierSpline, θ::Real, t::Real)
    return P(θ)(t)
end

function (P::FourierSpline)(θ::Real, t::Real)
    return P(θ)(t)
end