module FourierSplines
export Spline, Knots, knots, coefficients, degree, order, dimension, eval, FourierSpline, truncate, norm

import LinearAlgebra

include("Knots.jl")
include("Splines.jl")
include("FourierSplinesClass.jl")
include("LinAlgExt.jl")

end