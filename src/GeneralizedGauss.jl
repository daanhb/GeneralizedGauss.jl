module GeneralizedGauss

using BasisFunctions, LinearAlgebra
import BasisFunctions: moment

export compute_moments,
    compute_gauss_rule, compute_gauss_rules,
    oscillatory_basis_exp, oscillatory_basis_cos, log_basis,
    jacobi_matrix, newton

import Base:
    eltype,
    length,
    size



include("rootfinding.jl")

include("sets.jl")

include("quadrule.jl")

include("representations.jl")

include("gengauss.jl")

end # module
