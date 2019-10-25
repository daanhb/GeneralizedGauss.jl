
compute_moments(basis::Dictionary) = [moment(basis,i) for i in eachindex(basis)]

function apply_quad(w, x, f)
    z = w[1] * f(x[1])
    for i in 2:length(w)
        z += w[i] * f(x[i])
    end
    z
end

compute_weights(x, basis, B) = interpolation_matrix(basis, x)' \ B

"""
A container to collect information about a quadrature rule, mainly for use
in types describing exactness equations.
"""
struct QuadRuleData{B <: Dictionary,T}
    "The basis for which the rule is exact."
    basis   ::  B
    "The length of the quadrature rule"
    len     ::  Int
    "Moments of the basis with respect to some measure."
    moments ::  Vector{T}
end

eltype(::QuadRuleData{B,T}) where {B,T} = T

function QuadRuleData(basis::Dictionary{S,T}, len::Int, moments = compute_moments(basis)) where {S,T}
    QuadRuleData{typeof(basis),T}(basis, len, moments)
end

length(d::QuadRuleData) = d.len

basis(d::QuadRuleData) = d.basis

moments(d::QuadRuleData) = d.moments



"A type that describes the nonlinear equations of exactness of a quadrature rule, for use in Newton's method."
abstract type QuadRuleSystem{T} <: NonlinearSystem end

eltype(::QuadRuleSystem{T}) where {T} = T

data(sys::QuadRuleSystem) = sys.data

for op in [:basis, :moments]
    @eval $op(sys::QuadRuleSystem) = $op(sys.data)
end

length(rule::QuadRuleSystem) = dofs(rule)

size(rule::QuadRuleSystem) = (length(rule), dofs(rule))

quadlength(rule::QuadRuleSystem) = length(rule.data)


# Evaluate the residual of the system of equations
function residual!(result, sys::QuadRuleSystem, newton_x)
    w,x = newton_to_quad(sys, newton_x)
    residual!(result, sys, w, x, basis(sys), moments(sys))
end

function residual!(result, sys::QuadRuleSystem, w, x, basis, moments)
    for (i,ϕ) in enumerate(basis)
        result[i] =  apply_quad(w, x, ϕ) - moments[i]
    end
    result
end

# Evaluate the Jacobian of the system of equations
function jacobian!(J, sys::QuadRuleSystem, newton_x)
    w,x = newton_to_quad(sys, newton_x)
    jacobian!(J, sys, w, x, basis(sys))
end


function newton_to_quad(sys::QuadRuleSystem, newton_x)
    T = eltype(newton_x)
    x = zeros(T, length(sys.data))
    w = zeros(T, length(sys.data))
    newton_to_quad!(sys, w, x, newton_x)
end

function quad_to_newton(sys::QuadRuleSystem, w, x)
    T = promote_type(eltype(w),eltype(x))
    newton_x = zeros(T, dofs(sys))
    quad_to_newton!(sys, w, x, newton_x)
end


"System of equations for a quadrature rule where all points are free."
struct QuadRuleFreePoints{B,T} <: QuadRuleSystem{T}
    data    ::  QuadRuleData{B,T}
end

"The total number of degrees of freedom in the system."
dofs(rule::QuadRuleFreePoints) = 2*length(rule.data)


"System of equations for a quadrature rule with some points fixed."
struct QuadRuleFixedPoints{B,T} <: QuadRuleSystem{T}
    "Data about the quadrature rule: basis, moments, etcetera."
    data        ::  QuadRuleData{B,T}
    "Indices of the fixed points."
    fixed_idxs  ::  Vector{Int}
    "Values of the fixed points."
    fixed_pts   ::  Vector{T}
    "Indices of the free points."
    free_idxs   ::  Vector{Int}

    function QuadRuleFixedPoints{B,T}(data, fixed_idxs, fixed_pts) where {B,T}
        P = length(fixed_idxs)
        @assert length(fixed_pts) == P
        @assert P > 0
        # Make sure indices and points are sorted.
        for i in 1:P-1
            @assert fixed_idxs[i] < fixed_idxs[i+1]
            @assert fixed_pts[i] < fixed_pts[i+1]
        end
        new(data, fixed_idxs, fixed_pts, setdiff(1:length(data), fixed_idxs))
    end
end

QuadRuleFixedPoints(data::QuadRuleData{B,T}, fixed_idxs, fixed_pts) where {B,T} =
    QuadRuleFixedPoints{B,T}(data, fixed_idxs, fixed_pts)

"The number of fixed points in the quadrature rule."
nbfixed(rule::QuadRuleFixedPoints) = length(rule.fixed_pts)

"The number of free points in the quadrature rule."
nbfree(rule::QuadRuleFixedPoints) = length(rule.free_idxs)

dofs(rule::QuadRuleFixedPoints) = 2*length(rule.data) - nbfixed(rule)

function setfixedpoint!(rule::QuadRuleFixedPoints, idx, xstar)
    rule.fixed_pts[rule.fixed_idxs[idx]] = xstar
end

# newton_x contains l weights, followed by l points
function newton_to_quad!(sys::QuadRuleFreePoints, w, x, newton_x)
    @assert length(newton_x) == dofs(sys)
    l = length(w)
    for i in 1:l
        w[i] = newton_x[i]
        x[i] = newton_x[l+i]
    end
    w,x
end

function quad_to_newton!(sys::QuadRuleFreePoints, w, x, newton_x)
    @assert length(newton_x) == dofs(sys)
    l = length(w)
    for i in 1:l
        newton_x[i] = w[i]
        newton_x[l+i] = x[i]
    end
    newton_x
end


# newton_x contains l weights, followed by l-nbfixed(rule) free points (the other points are fixed)
function newton_to_quad!(sys::QuadRuleFixedPoints, w, x, newton_x)
    @assert length(newton_x) == dofs(sys)

    l = length(w)
    P = nbfixed(sys)

    for i in 1:l
        w[i] = newton_x[i]
    end

    # Copy the free points into x
    for (j,i) in enumerate(sys.free_idxs)
        x[i] = newton_x[l+j]
    end
    # and copy the fixed points into x
    for (i,xi) in zip(sys.fixed_idxs, sys.fixed_pts)
        x[i] = xi
    end
    w,x
end

function quad_to_newton!(sys::QuadRuleFixedPoints, w, x, newton_x)
    @assert length(newton_x) == dofs(sys)

    l = length(w)
    P = nbfixed(sys)

    for i in 1:l
        newton_x[i] = w[i]
    end

    # Copy the free points
    for (j,i) in enumerate(sys.free_idxs)
        newton_x[l+j] = x[i]
    end
    newton_x
end

function jacobian!(J, sys::QuadRuleFreePoints, w, x, basis)
    l = length(w)

    # We are computing the derivative of \sum_{i=1}^l w_i ϕ_j(x_i) wrt w_i and x_i
    # - First the derivatives with respect to the weight w[i]
    for j in 1:length(basis)
        for i in 1:l
            J[j,i] = eval_element(basis, j, x[i])
        end
    end
    # - Then the weight times the derivatives of the basis functions
    for j in 1:length(basis)
        for i in 1:l
            J[j,l+i] = w[i] * eval_element_derivative(basis, j, x[i])
        end
    end
    J
end


function jacobian!(J, sys::QuadRuleFixedPoints, w, x, basis)
    l = length(w)
    P = nbfixed(sys)

    # We are computing the derivative of \sum_{i=1}^l w_i ϕ_j(x_i) wrt w_i and x_i
    # as before, but now some of the x_i's are fixed and we have to skip them.
    # - First the derivatives with respect to the weight w[i]: same as before
    for j in 1:length(basis)
        for i in 1:l
            J[j,i] = eval_element(basis, j, x[i])
        end
    end
    # - Then the weight times the derivatives of the basis functions
    for j in 1:length(basis)
        for (k,i) in enumerate(sys.free_idxs)
            J[j,l+k] = w[i] * eval_element_derivative(basis, j, x[i])
        end
    end
    J
end
