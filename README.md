# GeneralizedGauss.jl

A package for the computation of generalized Gaussian quadrature rules. For a description of the algorithm, see the paper [On the computation of Gaussian quadrature rules for Chebyshev sets of linearly independent functions](https://arxiv.org/abs/1710.11244).

First of all, the package can be used to compute classical Gaussian quadrature rules. This is not recommended: the standard methods using orthogonal polynomials and the Jacobi matrix eigenvalue problem is much more efficient! But we can use it to illustrate the algorithm:
```julia
julia> using BasisFunctions, DomainSets, GeneralizedGauss

julia> basis = ChebyshevT(10)
ChebyshevT(10)

julia> w, x = compute_gauss_rule(basis)
([0.23692688505618958, 0.47862867049936764, 0.5688888888888858, 0.4786286704993679, 0.23692688505618886], [-0.9061798459386639, -0.5384693101056822, 1.3425672027729429e-15, 0.538469310105682, 0.9061798459386639])
```
This is the classical Gauss-Legendre rule with 5 points on `[-1,1]`. The quadrature rule depends on the space, not on the chosen basis for that space, hence we could have used any basis of polynomials. Here, we used Chebyshev polynomials. Of course these 5 quadrature points are the roots of the Legendre polynomial of degree 5.

The method computes a sequence of quadrature rules of increasing length. All intermediate rules can be obtained by invoking `compute_gauss_rules` instead of `compute_gauss_rule`. Also, informative output is displayed by passing the `verbose` flag.

```julia
julia> w, x, xi_upper, xi_lower, w_lower, w_upper, x_lower, x_upper = compute_gauss_rules(basis, verbose=true);
Computing initial one point rule
One point quadrature rule is: [0.0], [2.0]
Estimating upper canonical representation near 0.0
Upper principal representation 1 : xi is -0.3333333333333333
Estimating lower canonical representation near -0.3333333333333333
Lower principal representation 2 : xi is -0.5773502691896282
Estimating upper canonical representation near -0.5773502691896282
Upper principal representation 2 : xi is -0.6898979485563844
Estimating lower canonical representation near -0.6898979485563844
Lower principal representation 3 : xi is -0.774596669239068
Estimating upper canonical representation near -0.774596669239068
Upper principal representation 3 : xi is -0.8228240809745913
Estimating lower canonical representation near -0.8228240809745913
Lower principal representation 4 : xi is -0.8611363115939351
Estimating upper canonical representation near -0.8611363115939351
Upper principal representation 4 : xi is -0.8857916077709647
Estimating lower canonical representation near -0.8857916077709647
Lower principal representation 5 : xi is -0.9061798459386639

julia> x_lower
5-element Vector{Vector{Float64}}:
 [0.0]
 [-0.5773502691896282, 0.5773502691896679]
 [-0.774596669239068, 1.6279410648108757e-11, 0.7745966692339527]
 [-0.8611363115939351, -0.33998104358429126, 0.33998104358465014, 0.8611363115940404]
 [-0.9061798459386639, -0.5384693101056822, 1.3425672027729429e-15, 0.538469310105682, 0.9061798459386639]

julia> x_upper
4-element Vector{Vector{Float64}}:
 [-0.3333333333333333, 1.0]
 [-0.6898979485563844, 0.2898979485568163, 1.0]
 [-0.8228240809745913, -0.1810662711185286, 0.5753189235216944, 1.0]
 [-0.8857916077709647, -0.44631397272375223, 0.1671808647378336, 0.720480271312439, 1.0]
```
The method returns upper and lower principal representations. It is describe in the above-mentioned paper what they are. In this case, the lower principal representations in `x_lower` are the Gauss-Legendre quadrature points for `n=1,2,3,4,5`. The upper principal representations correspond to Gauss-Radau rules, which include the right endpoint `x=1`.

Quadrature rules can be computed for more general function spaces as well. For this we use the [BasisFunctions](https://github.com/JuliaApproximation/BasisFunctions.jl) package.
```julia
using BasisFunctions, DomainSets, GeneralizedGauss
N = 6
cheb = ChebyshevT(N) → 0..1        # Chebyshev polynomials scaled to [0,1]
basis = cheb ⊕ log*cheb                # space of the form p(x) + log(x)*q(x) with p and q polynomials
w, x = compute_gauss_rule(basis)
```
The result of this computation is a quadrature rule with 6 points, that will integrate polynomials up to degree 5 exactly on `[0,1]`, as well as `log(x)` times the same set of polynomials.

The rule can be used to integrate functions with a logarithmic singularity of the form `f(x) = p(x) + log(x) q(x)`, without having to know what `p` and `q` are since only `f` itself is evaluated.

One example is the Hankel function of the first kind and order zero, which has a logarithmic singularity at `x=0`. We can compare against an adaptive numerical integration procedure such as the one provide by the QuadGK package:
```julia
julia> using QuadGK, SpecialFunctions

julia> quadgk(x->besselh(0, 1, x), 0, 1)[1]
0.9197304100897603 - 0.6370693745520851im

julia> sum(w[k]*besselh(0, 1, x[k]) for k in 1:6)
0.919730410083813 - 0.6370693754207613im
```

The method in this package is not terribly efficient. For alternative computational schemes, consider the papers on this topic by V. Rokhlin, J. Bremer and others on generalized Gaussian quadrature.
