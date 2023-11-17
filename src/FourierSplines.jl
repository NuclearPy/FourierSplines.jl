module FourierSplines
import LinearAlgebra

include("Knots.jl")
include("Splines.jl")
export Spline, Knots, knots, coefficients, degree, order, dimension, evaluate, FourierSpline, truncate, indices

include("FourierSplinesClass.jl")
export FourierSpline, truncate

include("LinAlgExt.jl")
export SplineLinOp, norm, inv, det

end