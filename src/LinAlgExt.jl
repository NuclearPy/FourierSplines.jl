# Norms
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