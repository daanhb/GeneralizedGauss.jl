# GeneralizedGauss.jl

A package for the computation of generalized Gaussian quadrature rules. For a description of the algorithm, see the paper [On the computation of Gaussian quadrature rules for Chebyshev sets of linearly independent functions](https://arxiv.org/abs/1710.11244).

## Classical polynomial Gaussian quadrature

First of all, the package can be used to compute the standard Gauss-Legendre quadrature rules. This is not recommended, as existing methods are much more efficient. Yet, we can use it to illustrate the syntax and the algorithm:
```julia
julia> using BasisFunctions, DomainSets, GeneralizedGauss

julia> basis = ChebyshevT(10)
ChebyshevT(10)

julia> w, x = compute_gauss_rule(basis)
([0.23692688505618958, 0.47862867049936764, 0.5688888888888858, 0.4786286704993679, 0.23692688505618886], [-0.9061798459386639, -0.5384693101056822, 1.3425672027729429e-15, 0.538469310105682, 0.9061798459386639])
```
The result is the classical Gauss-Legendre rule with 5 points on `[-1,1]`. That quadrature rule depends on the space of polynomials, and not on the chosen basis. Here, we have arbitrarily used a basis of Chebyshev polynomials. Of course the resulting 5 quadrature points are the roots of the Legendre polynomial of degree 5.

The method of this package computes a sequence of quadrature rules of increasing length. All intermediate rules can be obtained by invoking `compute_gauss_rules` instead of `compute_gauss_rule`. Also, informative output is displayed by passing the `verbose` flag.

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
The method returns upper and lower principal representations. It is described in the above-mentioned paper what they are. In this case, the lower principal representations in `x_lower` are the Gauss-Legendre quadrature points for `n=1,2,3,4,5`. The upper principal representations correspond to Gauss-Radau rules, which include the right endpoint `x=1`.

## A function space with a logarithmic singularity

Quadrature rules can be computed for more general function spaces as well. For this we use the [BasisFunctions](https://github.com/JuliaApproximation/BasisFunctions.jl) package.
```julia
using BasisFunctions, DomainSets, GeneralizedGauss
N = 6
cheb = ChebyshevT(N) → 0..1        # Chebyshev polynomials scaled to [0,1]
basis = cheb ⊕ log*cheb            # space of the form p(x) + log(x)*q(x) with p and q polynomials
w, x = compute_gauss_rule(basis)
```
The result of this computation is a quadrature rule with 6 points, that will integrate polynomials up to degree 5 exactly on `[0,1]`, as well as `log(x)` times the same set of polynomials. The rule can be used to integrate functions with a logarithmic singularity of the form `f(x) = p(x) + log(x) q(x)`, without having to know what `p` and `q` are since only `f` itself is evaluated.

One example is the Hankel function of the first kind and order zero, which has a logarithmic singularity at `x=0`. We can compare against an adaptive numerical integration procedure such as the one provide by the [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package:
```julia
julia> using QuadGK, SpecialFunctions

julia> quadgk(x->besselh(0, 1, x), 0, 1)[1]
0.9197304100897603 - 0.6370693745520851im

julia> sum(w[k]*besselh(0, 1, x[k]) for k in 1:6)
0.919730410083813 - 0.6370693754207613im
```

The basis of this example rapidly becomes ill-conditioned when the polynomial degree is increased. Computations should be performed in higher-precision arithmetic, or a basis with better conditioning should be used -- for example by orthogonalizing in the space.

## Defining your own basis functions

The Gauss-Legendre rule above can also be computed by supplying polynomials and their derivatives "by hand":
```julia
n = 10
funs = [x->x^i for i in 0:n-1]
fun_derivs = vcat(x->zero(x), [x->i*x^(i-1) for i in 1:n-1])
basis = quadbasis(funs, fun_derivs, -1.0, 1.0)
w, x = compute_gauss_rule(basis)

([0.23692688618936705, 0.4786286821172505, 0.5688889121651178, 0.4786286156288994, 0.23692690389936527], [-0.9061798450949421, -0.5384693038382261, 2.5128577433358772e-8, 0.5384693147514663, 0.9061798400098149])
```
Or, in high precision:
```julia
T = BigFloat
basis = quadbasis(funs, fun_derivs, -one(T), one(T))
w, x = compute_gauss_rule(basis)

(BigFloat[0.2369268850561890875142640407199173626432600022124140155828612927606881455227749, 0.4786286704993664680412915148356381929122955533431415399727531757892658464848098, 0.5688888888888888888888888888888888888888888888888888888887562373626443396229523, 0.4786286704993664680412915148356381929122955533431415399727437693492163414736387, 0.2369268850561890875142640407199173626432600022124140155828855247381853268958114], BigFloat[-0.9061798459386639927976268782993929651256519107625308628737420458809403734172487, -0.538469310105683091036314420700208804967286606905559956202220239640514388749385, 9.586940463475928559453206027487843641529415939603793230907119240680334648035902e-60, 0.5384693101056830910363144207002088049672866069055599562021960699000907133035065, 0.9061798459386639927976268782993929651256519107625308628737391095560040519914682])```

## Other references

The method in this package is not at all efficient. For alternative computational schemes, consider the papers on this topic by V. Rokhlin, J. Bremer and others on generalized Gaussian quadrature.
