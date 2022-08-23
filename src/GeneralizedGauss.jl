module GeneralizedGauss

using BasisFunctions, LinearAlgebra, NLsolve

import BasisFunctions: moment

export quadbasis,
    compute_moments,
    compute_gauss_rule,
    compute_gauss_rules

import Base:
    eltype,
    length,
    size

include("basis.jl")
include("quadrule.jl")
include("representations.jl")
include("gengauss.jl")

end # module
