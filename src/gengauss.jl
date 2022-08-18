
Fk(w, x, u2k, c2k) = apply_quad(w, x, u2k) - c2k
Gk(w, x, u2kp1, c2kp1) = apply_quad(w, x, u2kp1) - c2kp1

default_threshold(dict::Dictionary) = default_threshold(codomaintype(dict))
default_threshold(::Type{T}) where {T <: AbstractFloat} = sqrt(eps(T))/100
default_threshold(::Type{BigFloat}) = big(1e-20)

solver_tolerance(::Type{Float64}) = 1e-8
solver_tolerance(::Type{T}) where {T} = sqrt(eps(T))
solver_tolerance(::Type{BigFloat}) = BigFloat(1e-30)

function ismonotonic(values)
    if length(values) <= 1
        return true
    end
    if first(values) > last(values)
        reduce(&, values[1:end-1] .> values[2:end])
    else
        reduce(&, values[1:end-1] .< values[2:end])
    end
end

function solve_system(rule, w0, x0; verbose=false, options...)
    x_init = quad_to_newton(rule, w0, x0)
    F!(Fx, x) = residual!(Fx, rule, x)
    J!(Jx, x) = jacobian!(Jx, rule, x)

    tol = solver_tolerance(eltype(x_init))
    r = nlsolve(F!, J!, x_init; ftol = tol, options...)
    w, x = newton_to_quad(rule, r.zero)
    converged(r), w, x
end

function supportleft(dict)
    t = leftendpoint(support(dict))
    isinf(t) ? -one(t)*10 : t
end
supportright(dict) = rightendpoint(support(dict))

function compute_one_point_rule(dict, moments; options...)
    @assert length(dict) == 2
    @assert length(moments) == 2

    x0 = 1/2 * (supportleft(dict) + supportright(dict))
    w0 = moments[1] / eval_element(dict, 1, x0)
    converged, w, x = compute_lower_principal_representation(dict, moments, w0, x0; options...)
    @assert converged
    w, x
end


function compute_upper_canonical_representation(dict, moments, xi, w0, x0;
        verbose=false, options...)
    @assert iseven(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)>>1)+1

    rule = CanonicalRepresentationOdd_K1(dict, xi, moments)
    try
        solve_system(rule, w0, x0; verbose, options...)
    catch e
        if e isa InterruptException
            rethrow()
        end
        println("ERROR THROWN at $(xi) in computation of upper canonical")
        @show e
        false, w0, x0
    end
end


function compute_many_upper_canonical_representation(dict, moments, w0, x0, pts;
        verbose=false, options...)

    l = length(dict) >> 1
    n = length(pts)
    T = eltype(w0)
    w = zeros(T,l+1,n)
    x = zeros(T,l+1,n)
    hasconverged = zeros(Bool,n)
    for i = n:-1:1
        xi = pts[i]
        if i < n
            w0 = w[:,i+1]
            x0 = x[:,i+1]
        end
        converged, w1, x1 = compute_upper_canonical_representation(dict, moments, xi, w0, x0; verbose, options...)
        if !converged && verbose
            println("Many upper canonical: not converged for $(xi)")
        end
        w[:,i] = w1
        x[:,i] = x1
        hasconverged[i] = converged
    end
    I = findall(hasconverged)
    any(hasconverged), w[:,I], x[:,I], pts[I]
end

function estimate_upper_canonical_representation(dict, moments, a, b, w0, x0;
        verbose=false, options...)
    @assert isodd(length(dict))
    @assert length(dict) == length(moments)
    l = (length(dict)-1) >> 1
    @assert length(w0) == l+1
    @assert length(x0) == l+1

    verbose && println("Estimating upper canonical representation near $(b)")
    n = 8
    estimate_upper_canonical_representation(dict, moments, a, b, w0, x0, n;
        verbose, options...)
end

function interpolate_starting_values(a, b, p, w_left, x_left, w_right, x_right)
    θ = (p-a)/(b-a)
    w0 = w_left + θ * (w_right-w_left)
    x0 = x_left + θ * (x_right-x_left)
    w0, x0
end

switches_sign(values) = ! (all(values .> 0) || all(values .< 0))

function estimate_upper_canonical_representation(dict, moments, a, b, wb, xb, n, wa=nothing, xa=nothing;
        verbose=false, options...)
    @assert n <= 1024

    l = (length(dict)-1) >> 1
    pts = range(a, b, length=n+2)[2:end-1]
    if !(wa == nothing)
        # w0, x0 = interpolate_starting_values(a, b, pts[end], wa, xa, wb, xb)
        w0, x0 = wa, xa
    else
        w0, x0 = wb, xb
    end

    some_converged, w, x, pts2 = compute_many_upper_canonical_representation(dict[1:2*l], moments[1:2*l], w0, x0, pts; verbose, options...)
    if some_converged
        if verbose && length(pts2) < length(pts)
            println("Upper canonical: converged for $(length(pts2)) out of $(length(pts)) points")
        end
        Fvals = [Fk(w[:,i],x[:,i],dict[2*l+1], moments[2*l+1]) for i in 1:size(w,2)]
        if verbose && !ismonotonic(Fvals)
            println("Upper canonical: function F is not monotonic but should be.")
        end
        if switches_sign(Fvals)
            if first(Fvals) > 0
                I = findlast(Fvals .> 0)
            else
                I = findlast(Fvals .< 0)
            end
            pts2[I], pts2[I+1], w[:,I], x[:,I], w[:,I+1], x[:,I+1]
        else
            if first(Fvals) > 0
                if first(Fvals) > last(Fvals)
                    a_new = pts2[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = w0; xb_new = x0;
                else
                    a_new = a; b_new = pts2[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                end
            else
                if first(Fvals) > last(Fvals)
                    a_new = a; b_new = pts2[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                else
                    a_new = pts2[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = w0; xb_new = x0;
                end
            end
            verbose && println("Upper canonical: refining from $((a,b)) to $((a_new,b_new))")
            estimate_upper_canonical_representation(dict, moments, a_new, b_new, wb_new, xb_new, n, wa_new, xa_new; verbose, options...)
        end
    else
        verbose && println("Upper canonical: increasing n to $(2n)")
        estimate_upper_canonical_representation(dict, moments, a, b, wb, xb, 2n; verbose, options...)
    end
end

function compute_upper_principal_representation(dict, moments, w0, x0;
        verbose=false, options...)
    @assert isodd(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)>>1)+1

    rule = UpperPrincipalEven(dict, moments)
    try
        solve_system(rule, w0, x0; verbose, options...)
    catch e
        if e isa InterruptException
            rethrow()
        end
        println("ERROR THROWN in computation of upper principal")
        @show e
        false, w0, x0
    end
end

function compute_lower_canonical_representation(dict, moments, xi, w0, x0;
        verbose=false, options...)
    @assert isodd(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)+1)>>1

    rule = CanonicalRepresentationEven_J1(dict, xi, moments)
    try
        solve_system(rule, w0, x0; verbose, options...)
    catch e
        if e isa InterruptException
            rethrow()
        end
        println("ERROR THROWN at $(xi) in computation of lower canonical")
        @show e
        false, w0, x0
    end
end

function compute_many_lower_canonical_representation(dict, moments, w0, x0, pts;
        verbose=false, options...)

    l = (length(dict)+1) >> 1
    n = length(pts)
    T = eltype(w0)
    w = zeros(T,l,n)
    x = zeros(T,l,n)
    hasconverged = zeros(Bool,n)
    for i = n:-1:1
        xi = pts[i]
        if i < n
            w0 = w[:,i+1]
            x0 = x[:,i+1]
        end
        converged, w1, x1 =
            compute_lower_canonical_representation(dict, moments, xi, w0, x0;
                verbose, options...)
        if !converged && verbose
            println("Many lower canonical: not converged for $(xi)")
        end
        w[:,i] = w1
        x[:,i] = x1
        hasconverged[i] = converged
    end
    I = findall(hasconverged)
    any(hasconverged), w[:,I], x[:,I], pts[I]
end

function estimate_lower_canonical_representation(dict, moments, a, b, w0, x0;
        verbose = false, options...)
    @assert iseven(length(dict))
    @assert length(dict) == length(moments)
    l = length(dict) >> 1
    @assert length(w0) == l
    @assert length(x0) == l

    verbose && println("Estimating lower canonical representation near $(b)")
    n = 8
    estimate_lower_canonical_representation(dict, moments, a, b, w0, x0, n; verbose, options...)
end

function estimate_lower_canonical_representation(dict, moments, a, b, wb, xb, n, wa=nothing, xa=nothing;
        verbose=false, options...)
    @assert n <= 1024
    l = length(dict) >> 1
    pts = range(a, b, length=n+2)[2:end-1]
    if !(wa==nothing)
        w0, x0 = wa, xa
    else
        w0, x0 = wb, xb
    end
    someconverged, w, x, pts2 = compute_many_lower_canonical_representation(dict[1:2*l-1], moments[1:2*l-1], w0, x0, pts; verbose, options...)
    if someconverged
        if verbose && length(pts2)<length(pts)
            println("Lower canonical: converged for $(length(pts2)) out of $(length(pts)) points")
        end
        Fvals = [Fk(w[:,i],x[:,i],dict[2*l], moments[2*l]) for i in 1:size(w,2)]
        if verbose && !ismonotonic(Fvals)
            println("Lower canonical: function F is not monotonic but should be.")
        end
        if (first(Fvals)*last(Fvals) < 0)
            if first(Fvals) > 0
                I = findlast(Fvals .> 0)
            else
                I = findlast(Fvals .< 0)
            end
            pts2[I], pts2[I+1], w[:,I], x[:,I], w[:,I+1], x[:,I+1]
        else
            if first(Fvals) > 0
                if first(Fvals) > last(Fvals)
                    a_new = pts2[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = wb; xb_new = xb;
                else
                    a_new = a; b_new = pts2[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                end
            else
                if first(Fvals) > last(Fvals)
                    a_new = a; b_new = pts2[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                else
                    a_new = pts2[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = wb; xb_new = xb;
                end
            end
            verbose && println("Lower canonical: refining from $((a,b)) to $((a_new,b_new))")
            estimate_lower_canonical_representation(dict, moments, a_new, b_new, wb_new, xb_new, 8, wa_new, xa_new; verbose, options...)
        end
    else
        verbose && println("Lower canonical: no convergence on grid, increasing to n=$(2n)")
        estimate_lower_canonical_representation(dict, moments, a, b, wb, xb, 2n, wa, xa; verbose, options...)
    end
end

function compute_lower_principal_representation(dict, moments, w0, x0;
        verbose=false, options...)
    rule = LowerPrincipalOdd(dict, moments)
    try
        solve_system(rule, w0, x0; verbose, options...)
    catch e
        if e isa InterruptException
            rethrow()
        end
        println("ERROR THROWN in computation of lower principal")
        @show e
        false, w0, x0
    end
end


function compute_gauss_rules(dict::Dictionary, moments = compute_moments(dict);
        verbose = false, options...)
    @assert iseven(length(dict))

    n = length(dict)
    l = n >> 1
    T = codomaintype(dict)

    xk_upper = fill(T[],0)
    wk_upper = fill(T[],0)
    xk_lower = fill(T[],0)
    wk_lower = fill(T[],0)

    verbose && println("Computing initial one point rule")
    wk,xk = compute_one_point_rule(dict[1:2], moments[1:2]; verbose, options...)
    verbose && println("One point quadrature rule is: ", xk, ", ", wk)

    xi_lower = zeros(T,l)
    xi_upper = zeros(T,l-1)
    xi_lower[1] = xk[1]

    push!(wk_lower, wk)
    push!(xk_lower, xk)

    if l == 1
        return wk, xk, xi_upper, xi_lower, wk_lower, wk_upper, xk_lower, xk_upper
    end

    for k = 1:l-1
        w0 = [wk; zero(T)]
        x0 = [xk; supportright(dict)]
        a = supportleft(dict)
        b = xk[1]
        p1, p2, w1, x1, w2, x2 =
            estimate_upper_canonical_representation(dict[1:2*k+1], moments[1:2*k+1], a, b, w0, x0; verbose, options...)
        converged, wk, xk = compute_upper_principal_representation(dict[1:2*k+1], moments[1:2*k+1], w2, x2; verbose, options...)
        @assert converged
        xi = xk[1]
        xi_upper[k] = xi
        verbose && println("Upper principal representation ", k, " : xi is ", xi)

        # Store the results for posterity
        push!(wk_upper, wk)
        push!(xk_upper, xk)

        w0 = wk
        x0 = xk
        a = supportleft(dict)
        b = xk[1]
        p1, p2, w1, x1, w2, x2 =
            estimate_lower_canonical_representation(dict[1:2*k+2], moments[1:2*k+2], a, b, w0, x0; verbose, options...)
        converged, wk, xk = compute_lower_principal_representation(dict[1:2*k+2], moments[1:2*k+2], w2, x2; verbose, options...)
        @assert converged
        xi = xk[1]
        xi_lower[k+1] = xi
        verbose && println("Lower principal representation ", k+1, " : xi is ", xi)

        push!(wk_lower, wk)
        push!(xk_lower, xk)
    end
    wk, xk, xi_upper, xi_lower, wk_lower, wk_upper, xk_lower, xk_upper
end

function compute_gauss_rule(dict::Dictionary, moments = compute_moments(dict); options...)
    wk, xk, xi_upper, xi_lower, wk_lower, wk_upper, xk_lower, xk_upper =
        compute_gauss_rules(dict, moments; options...)
    wk, xk
end
