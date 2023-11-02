LinearAlgebra.norm(S::Spline) = LinearAlgebra.norm(S.coefficients, Inf)
LinearAlgebra.norm(S::Spline, p::Real) = LinearAlgebra.norm(S.coefficients, p)