
Fk(w, x, u2k, c2k) = apply_quad(w, x, u2k) - c2k
Gk(w, x, u2kp1, c2kp1) = apply_quad(w, x, u2kp1) - c2kp1

default_threshold(dict::Dictionary) = default_threshold(codomaintype(dict))
default_threshold(::Type{T}) where {T <: AbstractFloat} = sqrt(eps(T))/100
default_threshold(::Type{BigFloat}) = big(1e-20)


supportleft(dict) = leftendpoint(support(dict))
supportright(dict) = rightendpoint(support(dict))

function compute_one_point_rule(dict, moments; verbose = false)
    @assert length(dict) == 2
    @assert length(moments) == 2

    x0 = 1/2 * (supportleft(dict) + supportright(dict))
    w0 = moments[1] / eval_element(dict, 1, x0)
    converged, w, x = compute_principal_representation_lower(dict, moments, w0, x0)
    @assert converged
    w, x
end


function compute_canonical_representation_upper(dict, moments, xi, w0, x0)
    @assert iseven(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)>>1)+1

    rule = CanonicalRepresentationOdd_K1(dict, xi, moments)
    x_init = quad_to_newton(rule, w0, x0)
    F!(Fx, x) = residual!(Fx, rule, x)
    J!(Jx, x) = jacobian!(Jx, rule, x)
    try
        r = nlsolve(F!, J!, x_init)
        w, x = newton_to_quad(rule, r.zero)
        converged(r), w, x
    catch e
        println("ERROR THROWN at $(xi) in computation of upper canonical")
        @show e
        false, w0, x0
    end
end


function compute_many_canonical_representation_upper(dict, moments, w0, x0, pts; verbose=false)
    local converged

    l = length(dict) >> 1
    n = length(pts)
    T = eltype(w0)
    w = zeros(T,l+1,n)
    x = zeros(T,l+1,n)
    for i = n:-1:1
        xi = pts[i]
        if i < n
            w0 = w[:,i+1]
            x0 = x[:,i+1]
        end
        converged, w1, x1 = compute_canonical_representation_upper(dict, moments, xi, w0, x0)
        if !converged
            verbose && println("Many upper canonical: not converged for $(xi)")
            break
        end
        w[:,i] = w1
        x[:,i] = x1
    end
    converged, w, x
end

function estimate_canonical_representation_upper(dict, moments, a, b, w0, x0; verbose=false)
    @assert isodd(length(dict))
    @assert length(dict) == length(moments)
    l = (length(dict)-1) >> 1
    @assert length(w0) == l+1
    @assert length(x0) == l+1

    verbose && println("Estimating upper canonical representation near $(b)")
    n = 8
    estimate_canonical_representation_upper(dict, moments, a, b, w0, x0, n; verbose)
end

function interpolate_starting_values(a, b, p, w_left, x_left, w_right, x_right)
    θ = (p-a)/(b-a)
    w0 = w_left + θ * (w_right-w_left)
    x0 = x_left + θ * (x_right-x_left)
    w0, x0
end

function estimate_canonical_representation_upper(dict, moments, a, b, wb, xb, n, wa=nothing, xa=nothing; verbose=false)
    @assert n <= 1024

    l = (length(dict)-1) >> 1
    pts = range(a, b, n+2)[2:end-1]
    if !(wa == nothing)
        w0, x0 = interpolate_starting_values(a, b, pts[end], wa, xa, wb, xb)
    else
        w0, x0 = wb, xb
    end

    converged, w, x = compute_many_canonical_representation_upper(dict[1:2*l], moments[1:2*l], w0, x0, pts; verbose)
    if converged
        Fvals = [Fk(w[:,i],x[:,i],dict[2*l+1], moments[2*l+1]) for i in 1:size(w,2)]
        ismonotonic(Fvals) || error("Upper canonical: function F is not monotonic but should be. Aborting.")
        if (first(Fvals)*last(Fvals) < 0)
            if first(Fvals) > 0
                I = findlast(Fvals .> 0)
            else
                I = findlast(Fvals .< 0)
            end
            pts[I], pts[I+1], w[:,I], x[:,I], w[:,I+1], x[:,I+1]
        else
            if first(Fvals) > 0
                if first(Fvals) > last(Fvals)
                    a_new = pts[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = w0; xb_new = x0;
                else
                    a_new = a; b_new = pts[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                end
            else
                if first(Fvals) > last(Fvals)
                    a_new = a; b_new = pts[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                else
                    a_new = pts[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = w0; xb_new = x0;
                end
            end
            verbose && println("Upper canonical: refining from $((a,b)) to $((a_new,b_new))")
            estimate_canonical_representation_upper(dict, moments, a_new, b_new, wb_new, xb_new, n, wa_new, xa_new; verbose)
        end
    else
        verbose && println("Upper canonical: increasing n to $(2n)")
        estimate_canonical_representation_upper(dict, moments, a, b, wb, xb, 2n; verbose)
    end
end

function compute_principal_representation_upper(dict, moments, w0, x0)
    @assert isodd(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)>>1)+1

    rule = UpperPrincipalEven(dict, moments)
    x_init = quad_to_newton(rule, w0, x0)
    F!(Fx, x) = residual!(Fx, rule, x)
    J!(Jx, x) = jacobian!(Jx, rule, x)
    try
        r = nlsolve(F!, J!, x_init)
        w, x = newton_to_quad(rule, r.zero)
        converged(r), w, x
    catch e
        println("ERROR THROWN in computation of upper principal")
        @show e
        false, w0, x0
    end
end

function compute_canonical_representation_lower(dict, moments, xi, w0, x0)
    @assert isodd(length(dict))
    @assert length(moments) == length(dict)
    @assert length(w0) == (length(dict)+1)>>1

    rule = CanonicalRepresentationEven_J1(dict, xi, moments)
    x_init = quad_to_newton(rule, w0, x0)
    F!(Fx, x) = residual!(Fx, rule, x)
    J!(Jx, x) = jacobian!(Jx, rule, x)
    try
        r = nlsolve(F!, J!, x_init)
        w, x = newton_to_quad(rule, r.zero)
        c = converged(r)
        converged(r), w, x
    catch e
        println("ERROR THROWN at $(xi) in computation of lower canonical")
        @show e
        false, w0, x0
    end
end

function compute_many_canonical_representation_lower(dict, moments, w0, x0, pts; verbose=false)
    local converged

    l = (length(dict)+1) >> 1
    n = length(pts)
    T = eltype(w0)
    w = zeros(T,l,n)
    x = zeros(T,l,n)
    for i = n:-1:1
        xi = pts[i]
        if i < n
            w0 = w[:,i+1]
            x0 = x[:,i+1]
        end
        converged, w1, x1 = compute_canonical_representation_lower(dict, moments, xi, w0, x0)
        if !converged
            verbose && println("Many lower canonical: not converged for $(xi)")
            break
        end
        w[:,i] = w1
        x[:,i] = x1
    end
    converged, w, x
end

function estimate_canonical_representation_lower(dict, moments, a, b, w0, x0; verbose = false)
    @assert iseven(length(dict))
    @assert length(dict) == length(moments)
    l = length(dict) >> 1
    @assert length(w0) == l
    @assert length(x0) == l

    verbose && println("Estimating lower canonical representation near $(b)")
    n = 8
    estimate_canonical_representation_lower(dict, moments, a, b, w0, x0, n; verbose)
end

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

function estimate_canonical_representation_lower(dict, moments, a, b, wb, xb, n, wa=nothing, xa=nothing; verbose=false)
    @assert n <= 1024
    l = length(dict) >> 1
    pts = range(a, b, n+2)[2:end-1]
    if !(wa==nothing)
        w0, x0 = interpolate_starting_values(a, b, pts[end], wa, xa, wb, xb)
    else
        w0 = wb
        x0 = xb
    end
    converged, w, x = compute_many_canonical_representation_lower(dict[1:2*l-1], moments[1:2*l-1], w0, x0, pts; verbose)
    if converged
        Fvals = [Fk(w[:,i],x[:,i],dict[2*l], moments[2*l]) for i in 1:size(w,2)]
        ismonotonic(Fvals) || error("Lower canonical: function F is not monotonic but should be. Aborting.")
        if (first(Fvals)*last(Fvals) < 0)
            if first(Fvals) > 0
                I = findlast(Fvals .> 0)
            else
                I = findlast(Fvals .< 0)
            end
            pts[I], pts[I+1], w[:,I], x[:,I], w[:,I+1], x[:,I+1]
        else
            if first(Fvals) > 0
                if first(Fvals) > last(Fvals)
                    a_new = pts[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = wb; xb_new = xb;
                else
                    a_new = a; b_new = pts[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                end
            else
                if first(Fvals) > last(Fvals)
                    a_new = a; b_new = pts[1];
                    wb_new = w[:,1]; xb_new = x[:,1];
                    wa_new = nothing; xa_new = nothing;
                else
                    a_new = pts[end]; b_new = b;
                    wa_new = w[:,end]; xa_new = x[:,end];
                    wb_new = wb; xb_new = xb;
                end
            end
            verbose && println("Lower canonical: refining from $((a,b)) to $((a_new,b_new))")
            estimate_canonical_representation_lower(dict, moments, a_new, b_new, wb_new, xb_new, n, wa_new, xa_new; verbose)
        end
    else
        verbose && println("Lower canonical: increasing n to $(2n)")
        estimate_canonical_representation_lower(dict, moments, a, b, wb, xb, 2n, wa, xa; verbose)
    end
end

function compute_principal_representation_lower(dict, moments, w0, x0)
    rule = LowerPrincipalOdd(dict, moments)
    x_init = quad_to_newton(rule, w0, x0)
    F!(Fx, x) = residual!(Fx, rule, x)
    J!(Jx, x) = jacobian!(Jx, rule, x)
    try
        r = nlsolve(F!, J!, x_init)
        w, x = newton_to_quad(rule, r.zero)
        converged(r), w, x
    catch e
        println("ERROR THROWN in computation of lower principal")
        @show e
        false, w0, x0
    end
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

    verbose && println("Computing initial one point rule")
    wk,xk = compute_one_point_rule(dict[1:2], moments[1:2]; verbose)
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
            estimate_canonical_representation_upper(dict[1:2*k+1], moments[1:2*k+1], a, b, w0, x0; verbose)
        converged, wk, xk = compute_principal_representation_upper(dict[1:2*k+1], moments[1:2*k+1], w2, x2)
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
            estimate_canonical_representation_lower(dict[1:2*k+2], moments[1:2*k+2], a, b, w0, x0; verbose)
        converged, wk, xk = compute_principal_representation_lower(dict[1:2*k+2], moments[1:2*k+2], w2, x2)
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
