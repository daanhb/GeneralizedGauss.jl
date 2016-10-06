module GeneralizedGauss

import Base: eltype, length, size

export oscillatory_basis_exp, oscillatory_basis_cos, log_basis, jacobi_matrix, newton

using BasisFunctions

import BasisFunctions: moment

include("rootfinding.jl")

include("sets.jl")

include("quadrule.jl")

include("representations.jl")

include("gengauss.jl")

end # module
