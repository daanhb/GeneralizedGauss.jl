
Fk(w, x, u2k, c2k) = apply_quad(w, x, u2k) - c2k

Gk(w, x, u2kp1, c2kp1) = apply_quad(w, x, u2kp1) - c2kp1

function leastsquares_approximate(dict::Dictionary, points, values)
    @assert length(points) == length(values)

    M = length(points)
    N = min(length(values), max(6,round(Int, sqrt(M))))
    basis = rescale(resize(dict, N), minimum(points), maximum(points))
    A = BasisFunctions.evaluation_matrix(basis, points)
    Expansion(basis, A\values)
end

default_threshold(dict::Dictionary) = default_threshold(codomaintype(dict))

default_threshold(::Type{T}) where {T <: AbstractFloat} = sqrt(eps(T))/100
default_threshold(::Type{BigFloat}) = big(1e-28)


supportleft(dict) = leftendpoint(support(dict))
supportright(dict) = rightendpoint(support(dict))

function compute_one_point_rule(dict, moments; verbose = false)
    @assert length(dict) == 2
    @assert length(moments) == 2

    x0 = 1/2 * (supportleft(dict) + supportright(dict))
    w0 = moments[1] / eval_element(dict, 1, x0)
    w, x = compute_principal_representation_lower(dict, moments, w0, x0)
    # w0 = supportright(dict) - supportleft(dict)
    # w, x, converged = compute_one_point_rule(dict, moments, x0, w0,
    #     maxiter=maxiter, verbose=verbose)
    # if !converged
    #     x0 = 1/2 * (supportleft(dict) + supportright(dict))
    #     w0 = supportright(dict) - supportleft(dict)
    # end
end

# function compute_one_point_rule(dict, moments, x0, w0; maxiter, verbose)
#     sys = LowerPrincipalOdd(dict, moments)
#     result = newton(sys, quad_to_newton(sys, w0, x0),
#         maxiter=maxiter, verbose=verbose)
#     iterations = result[2]
#     w, x = newton_to_quad(sys, result[1])
#     w, x, iterations<maxiter
# end


function compute_canonical_representation_upper(dict, moments, xi, w0, x0)
    @assert iseven(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)>>1)+1

    rule = CanonicalRepresentationOdd_K1(dict, xi, moments)
    nx, iter, normx = newton_with_restart(rule, quad_to_newton(rule, w0, x0),
        threshold = default_threshold(dict), verbose = false)
    w, x = newton_to_quad(rule, nx)
end


function compute_many_canonical_representation_upper(dict, moments, gw, gx, n,
    pts = linspace(supportleft(dict), gx[1], n))

    l = length(dict) >> 1
    w = zeros(l+1,n)
    x = zeros(l+1,n)
    for i = n:-1:1
        xi = pts[i]
        if i == n
            w0 = [gw; 0.0]
            x0 = [gx; supportright(dict)]
        else
            w0 = w[:,i+1]
            x0 = x[:,i+1]
        end
        w1,x1 = compute_canonical_representation_upper(dict, moments, xi, w0, x0)
        w[:,i] = w1
        x[:,i] = x1
    end
    w,x
end

function chebyshev_approximation(basis, vals)
    A = approximation_operator(basis)
    Expansion(basis, A*vals)
end

function approximate_canonical_representation_upper(dict, moments, gw, gx, n)
    @assert isodd(length(dict))
    @assert length(dict) == length(moments)
    l = (length(dict)-1) >> 1
    @assert length(gw) == l
    @assert length(gx) == l

    basis = ChebyshevT(n) → supportleft(dict)..gx[1]
    pts = interpolation_grid(basis)
    w,x = compute_many_canonical_representation_upper(dict[1:2*l], moments[1:2*l], gw, gx, n, pts)
    T = eltype(w)
    Fvals = T[Fk(w[:,i],x[:,i],dict[2*l+1], moments[2*l+1]) for i in 1:size(w,2)]
    F_fun = chebyshev_approximation(basis, Fvals)
    w_funs = [chebyshev_approximation(basis, w[i,:]') for i in 1:size(w,1)]
    x_funs = [chebyshev_approximation(basis, x[i,:]') for i in 1:size(x,1)]
    F_fun, w_funs, x_funs
end

function compute_principal_representation_upper(dict, moments, w0, x0)
    rule = UpperPrincipalEven(dict, moments)
    nx, iter, normx = newton_with_restart(rule, quad_to_newton(rule, w0, x0),
        threshold = default_threshold(dict), verbose = false)
    w, x = newton_to_quad(rule, nx)
end

function compute_canonical_representation_lower(dict, moments, xi, w0, x0)
    @assert isodd(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)+1)>>1

    rule = CanonicalRepresentationEven_J1(dict, xi, moments)
    nx, iter, normx = newton_with_restart(rule, quad_to_newton(rule, w0, x0),
        threshold = default_threshold(dict), verbose = false, maxstep = 10)
    w, x = newton_to_quad(rule, nx)
end

function compute_many_canonical_representation_lower(dict, moments, gw, gx, n,
    pts = linspace(supportleft(dict), gx[1], n))

    l = (length(dict)+1) >> 1
    w = zeros(l,n)
    x = zeros(l,n)
    w[:,n] = gw
    x[:,n] = gx
    for i = n:-1:1
        xi = pts[i]
        if i == n
            w0 = gw
            x0 = gx
        else
            w0 = w[:,i+1]
            x0 = x[:,i+1]
        end
        w1,x1 = compute_canonical_representation_lower(dict, moments, xi, w0, x0)
        w[:,i] = w1
        x[:,i] = x1
    end
    w,x
end

function approximate_canonical_representation_lower(dict, moments, gw, gx, n)
    @assert iseven(length(dict))
    @assert length(dict) == length(moments)
    l = length(dict) >> 1
    @assert length(gw) == l
    @assert length(gx) == l

    basis = ChebyshevT(n) → supportleft(dict)..gx[1]
    pts = interpolation_grid(basis)
    w,x = compute_many_canonical_representation_lower(dict[1:2*l-1], moments[1:2*l-1], gw, gx, n, pts)
    T = eltype(w)
    Fvals = T[Fk(w[:,i],x[:,i],dict[2*l], moments[2*l]) for i in 1:size(w,2)]
    F_fun = chebyshev_approximation(basis, Fvals)
    w_funs = [chebyshev_approximation(basis, w[i,:]') for i in 1:size(w,1)]
    x_funs = [chebyshev_approximation(basis, x[i,:]') for i in 1:size(x,1)]
    F_fun, w_funs, x_funs
end

function compute_principal_representation_lower(dict, moments, w0, x0)
    rule = LowerPrincipalOdd(dict, moments)
    nx, iter, normx = newton_with_restart(rule, quad_to_newton(rule, w0, x0),
        threshold = default_threshold(dict), verbose = false)
    w, x = newton_to_quad(rule, nx)
end


function compute_gauss_rules(dict::Dictionary, moments = compute_moments(dict); verbose = false)
    @assert iseven(length(dict))

    n = length(dict)
    l = n >> 1

    T = codomaintype(dict)

    xk_upper = fill(T[],0)
    wk_upper = fill(T[],0)
    xk_lower = fill(T[],0)
    wk_lower = fill(T[],0)
    f_funs_upper = Array{Any}(undef,0)
    w_funs_upper = Array{Any}(undef,0)
    x_funs_upper = Array{Any}(undef,0)
    f_funs_lower = Array{Any}(undef,0)
    w_funs_lower = Array{Any}(undef,0)
    x_funs_lower = Array{Any}(undef,0)

    verbose && println("Computing initial one point rule")
    wk,xk = compute_one_point_rule(dict[1:2], moments[1:2], verbose=verbose)
    verbose && println("One point quadrature rule is: ", xk, ", ", wk)
    xi_lower = zeros(T,l)
    xi_upper = zeros(T,l-1)
    xi_lower[1] = xk[1]

    push!(wk_lower, wk)
    push!(xk_lower, xk)

    if l == 1
        return wk, xk, xi_upper, xi_lower, wk_lower, wk_upper, xk_lower, xk_upper, f_funs_lower, f_funs_upper, w_funs_lower, w_funs_upper, x_funs_lower, x_funs_upper
    end

    bn = 16
    for k = 1:l-1
        f_fun, w_funs, x_funs = approximate_canonical_representation_upper(
            dict[1:2*k+1], moments[1:2*k+1], wk, xk, bn)
        xi0 = bisection(f_fun, supportleft(f_fun), supportright(f_fun); threshold = default_threshold(T))
        w0 = T[w_funs[k](xi0) for k in 1:length(x_funs)]
        x0 = T[x_funs[k](xi0) for k in 1:length(x_funs)]
        wk, xk = compute_principal_representation_upper(dict[1:2*k+1], moments[1:2*k+1], w0, x0)
        xi = xk[1]
        xi_upper[k] = xi
        verbose && println("Upper principal representation ", k, " : xi is ", xi)

        # Store the results for posterity
        push!(wk_upper, wk)
        push!(xk_upper, xk)
        push!(f_funs_upper, f_fun)
        push!(w_funs_upper, w_funs)
        push!(x_funs_upper, x_funs)

        f_fun, w_funs, x_funs = approximate_canonical_representation_lower(
            dict[1:2*k+2], moments[1:2*k+2], wk, xk, bn)
        xi0 = bisection(f_fun, supportleft(f_fun), supportright(f_fun); threshold = default_threshold(T))
        w0 = T[w_funs[k](xi0) for k in 1:length(w_funs)]
        x0 = T[x_funs[k](xi0) for k in 1:length(x_funs)]
        wk, xk = compute_principal_representation_lower(dict[1:2*k+2], moments[1:2*k+2], w0, x0)
        xi = xk[1]
        xi_lower[k+1] = xi
        if verbose
            println("Lower principal representation ", k+1, " : xi is ", xi)
        end

        push!(wk_lower, wk)
        push!(xk_lower, xk)
        push!(f_funs_lower, f_fun)
        push!(w_funs_lower, w_funs)
        push!(x_funs_lower, x_funs)

        residual = abs(sum([wk[i]*dict[2*k+2](xk[i]) for i=1:length(wk)]) - moments[2*k+2])

        if f_fun.coefficients[end] > max(default_threshold(T), 2*residual)
            verbose && println("Increasing bn parameter")
            bn *= 2
        end
    end
    wk, xk, xi_upper, xi_lower, wk_lower, wk_upper, xk_lower, xk_upper, f_funs_lower, f_funs_upper, w_funs_lower, w_funs_upper, x_funs_lower, x_funs_upper
end

function compute_gauss_rule(dict::Dictionary, moments = compute_moments(dict); options...)
    wk, xk, xi_upper, xi_lower, wk_lower, wk_upper, xk_lower, xk_upper, f_funs_lower, f_funs_upper, w_funs_lower, w_funs_upper, x_funs_lower, x_funs_upper =
        compute_gauss_rules(dict, moments; options...)
    wk, xk
end
