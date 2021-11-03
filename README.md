# GeneralizedGauss.jl

A package for the computation of generalized Gaussian quadrature rules. For a description of the algorithm, see the paper [On the computation of Gaussian quadrature rules for Chebyshev sets of linearly independent functions](https://arxiv.org/abs/1710.11244)

```julia
using BasisFunctions, DomainSets, GeneralizedGauss
N = 10
cheb = ChebyshevT(N) → 0..1        # Chebyshev polynomials scaled to [0,1]
basis = cheb ⊕ log*cheb                # space of the form p(x) + log(x)*q(x) with p and q polynomials
wk, xk = compute_gauss_rules(basis, verbose=true)
```

In this example, the majority of time is spent numerically computing the moments of the basis. Only afterwards the quadrature rule is being computed. In order to precompute the moments for this example, use:
```julia
function modified_moments(N, ::Type{T} = Float64) where {T}
    c = zeros(T, N+1)
    c[1] = 1
    for n in 2:(N>>1 + 1)
        c[2n-1] = -one(T) / ((2n-3)*(2n-1))
    end
    c
end
function modified_log_moments(N, ::Type{T} = Float64) where {T}
    chebA =  [BasisFunctions.rec_An(ChebyshevT{T}(N), k) for k in 0:N]
    chebB =  [BasisFunctions.rec_Bn(ChebyshevT{T}(N), k) for k in 0:N]
    chebC =  [BasisFunctions.rec_Cn(ChebyshevT{T}(N), k) for k in 0:N]
    modchebA = 2chebA
    modchebC = chebC
    modchebB = similar(chebB)
    for n in 1:length(chebB)
        modchebB[n] = chebB[n] - chebA[n]
    end
    L = zeros(T, N+1,N+1)
    L[1,:] = [logmoment(k,T) for k in 0:N]
    for n in 0:N-1
        for k in 0:N-1
            if n == 0
                L[n+2,k+1] = modchebA[n+1]*L[n+1,k+2] + modchebB[n+1]*L[n+1,k+1]
            else
                L[n+2,k+1] = modchebA[n+1]*L[n+1,k+2] + modchebB[n+1]*L[n+1,k+1] - modchebC[n+1]*L[n,k+1]
            end
        end
    end
    L[:,1]
end
```

The example becomes:
```julia
using BasisFunctions, DomainSets, GeneralizedGauss
N = 10
cheb = ChebyshevT(N) → 0..1
basis = cheb ⊕ log*cheb
moments = vcat(modified_moments(N-1), modified_log_moments(N-1))
wk, xk = compute_gauss_rules(basis, moments, verbose=true)
```
