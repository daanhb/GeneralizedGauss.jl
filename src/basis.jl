
"""
A generic set of functions.

This type represents a `Dictionary` as in the BasisFunctions package. The basis
functions are given by a vector of functions, and so are the derivatives.
The basis is defined on an interval `[a,b]`.
"""
struct GenericFunctionSet{S,T,F,D} <: Dictionary{S,T}
    funs        ::  F
    fun_derivs  ::  D
    a           ::  S
    b           ::  S

    function GenericFunctionSet{S,T,F,D}(funs::F, fun_derivs::D, a, b) where {S,T,F,D}
        if fun_derivs != nothing
            @assert length(funs) == length(fun_derivs)
        end
        new(funs, fun_derivs, a, b)
    end
end

GenericFunctionSet(funs, fun_derivs, a::S, b::T) where {S,T} =
    GenericFunctionSet(funs, fun_derivs, promote(a,b)...)
GenericFunctionSet(funs, fun_derivs, a::T, b::T) where {T} =
    GenericFunctionSet{T,T}(funs, fun_derivs, a, b)

GenericFunctionSet{S,T}(funs::F, fun_derivs::D, a, b) where {S,T,F,D} =
    GenericFunctionSet{S,T,F,D}(funs, fun_derivs, a, b)

Base.size(dict::GenericFunctionSet) = (length(dict.funs),)
BasisFunctions.support(dict::GenericFunctionSet) = dict.a..dict.b

BasisFunctions.unsafe_eval_element(basis::GenericFunctionSet, i, x) = basis.funs[i](x)
function BasisFunctions.unsafe_eval_element_derivative(basis::GenericFunctionSet, i, x, order)
    @assert order == 1
    basis.fun_derivs[i](x)
end


quadbasis(funs, fun_derivs, a, b) = GenericFunctionSet(funs, fun_derivs, a, b)


funeval(basis, i, x) = basis[i](x)
funeval(basis::Dictionary, i, x) =
    BasisFunctions.unsafe_eval_element(basis, i, x)
funeval_deriv(basis::Dictionary, i, x) =
    BasisFunctions.unsafe_eval_element_derivative(basis, i, x, 1)
