# Spline construction
mutable struct Spline{T <: Number}
    const knots :: Knots
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