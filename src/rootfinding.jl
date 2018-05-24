# rootfinding.jl

function bisection(fun, a, b, threshold = sqrt(eps(eltype(fun)))/100, maxiter = 100)
    root = (a+b)/2
    iter = 0
    while (abs(fun(root)) > threshold) && (iter < maxiter)
        if fun(a)*fun(root) < 0
            b = root
        else
            a = root
        end
        root = (a+b)/2
        iter += 1
    end
    if iter == maxiter
        println("Maximal number of iterations reached in bisection method.")
    end
    root
end



####################
# Newton iterations
#####################

"The abstract NonlinearSystem type represents a system of non-linear equations in n variables."
abstract type NonlinearSystem
end

# Default numeric type is Float64
# eltype(::Type{S}) where {S <: NonlinearSystem} = Float64

function residual(sys::NonlinearSystem, x)
    result = Array{eltype(sys)}(length(sys))
    residual!(result, sys, x)
end

function jacobian(sys::NonlinearSystem, x)
    J = Array{eltype(sys)}(size(sys))
    jacobian!(J, sys, x)
end



"""
Apply Newton iterations with the given starting values to find the roots of the given
system of equations. A convergence threshold and maximum number of iterations are optional.

The maximal step size in each iteration can be limited using `maxstep`.

Damped Newton iterations can be achieved by supplying damping parameter 0 < alpha < 1.
"""
function newton(sys::NonlinearSystem, x0; threshold = 1e-6, alpha = 1, maxiter = 100, maxstep = 1e5, verbose = false)
    n = length(x0)
    S = Array{eltype(sys)}(n)
    J = Array{eltype(sys)}(n, n)

    residual!(S, sys, x0)
    jacobian!(J, sys, x0)
    c = cond(J)
    threshold = max(threshold, eps(eltype(x0))*c)
    u = J \ S
    if norm(alpha*u) > maxstep
        u = u * maxstep/norm(u)
    end
    x1 = x0 - alpha*u

    iter = 1
    while norm(x1-x0) > threshold && iter < maxiter
        iter += 1
        x0 = x1
        residual!(S, sys, x0)
        jacobian!(J, sys, x0)
        u = J \ S
        if norm(alpha*u) > maxstep
            u = u * maxstep/norm(u)
        end
        x1 = x0 - alpha*u
    end
    if verbose
        if iter == maxiter
            println("Newton: maximum number of iterations, no convergence")
        else
            println("Newton: Converged to solution in ", iter, " iterations.")
        end
    end
    x1,iter,norm(x1-x0)
end


"Invoke Newton but restart with smaller maximal step size if an error is thrown."
function newton_with_restart(args...; maxstep = 1e5, maxiter = 100, options...)
    local nx,iter,normx
    nb_decreased = 0
    try
        nx,iter,normx = newton(args...; maxstep = maxstep, options...)
    catch e
        print("Error thrown by newton: ")
        println(e)
        iter = maxiter
    end
    while iter == maxiter
        println("Decreasing step size.")
        maxstep /= 2
        try
            nx,iter,normx = newton(args...; maxstep = maxstep, options...)
        catch e
            print("Error thrown by newton: ")
            println(e)
            iter = maxiter
        end
        nb_decreased += 1
        if nb_decreased > 20
            error("Too many step decreases in newton_with_restart. Aborting.")
        end
    end
    nx,iter,normx
end
