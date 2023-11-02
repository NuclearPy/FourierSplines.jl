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