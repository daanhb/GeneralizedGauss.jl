
"Generate the non-linear system for a lower principal representation for odd n."
function LowerPrincipalOdd(s::Dictionary, moments = compute_moments(s))
    n = length(s) - 1
    @assert isodd(n)
    l = (n+1) >> 1
    QuadRuleFreePoints(QuadRuleData(s, l, moments))
end

"Generate the non-linear system for an upper principal representation for odd n."
function UpperPrincipalOdd(s::Dictionary, moments = compute_moments(s))
    n = length(s) - 1
    @assert isodd(n)
    l = ((n+1) >> 1) + 1
    QuadRuleFixedPoints(QuadRuleData(s, l, moments), [1, l], [supportleft(s), supportright(s)])
end

"Generate the non-linear system for a lower principal representation for even n."
function LowerPrincipalEven(s::Dictionary, moments = compute_moments(s))
    n = length(s) - 1
    @assert iseven(n)
    l = (n >> 1) + 1
    QuadRuleFixedPoints(QuadRuleData(s, l, moments), [1], [supportleft(s)])
end

"Generate the non-linear system for an upper principal representation for even n."
function UpperPrincipalEven(s::Dictionary, moments = compute_moments(s))
    n = length(s) - 1
    @assert iseven(n)
    l = (n >> 1) + 1
    QuadRuleFixedPoints(QuadRuleData(s, l, moments), [l], [supportright(s)])
end

# The fixed root is in K1 = [a,t1].
function CanonicalRepresentationOdd_K1(s::Dictionary, xstar, moments = compute_moments(s))
    n = length(s) - 1
    @assert isodd(n)
    l = ((n+1) >> 1) + 1
    QuadRuleFixedPoints(QuadRuleData(s, l, moments), [1, l], [xstar, supportright(s)])
end

# The fixed root is in J1 = [t1,s2].
function CanonicalRepresentationOdd_J1(s::Dictionary, xstar, moments = compute_moments(s))
    n = length(s) - 1
    @assert isodd(n)
    l = ((n+1) >> 1) + 1
    QuadRuleFixedPoints(QuadRuleData(s, l, moments), [1, 2], [supportleft(s), xstar])
end

# The fixed root is in J1 = [a,s1].
function CanonicalRepresentationEven_J1(s::Dictionary, xstar, moments = compute_moments(s))
    n = length(s) - 1
    @assert iseven(n)
    l = (n >> 1) + 1
    QuadRuleFixedPoints(QuadRuleData(s, l, moments), [1], [xstar])
end

# The fixed root is in K1 = [s1,t2].
function CanonicalRepresentationEven_K1(s::Dictionary, xstar, moments = compute_moments(s))
    n = length(s) - 1
    @assert iseven(n)
    l = (n >> 1) + 2
    QuadRuleFixedPoints(QuadRuleData(s, l, moments), [1, 2, l], [supportleft(s), xstar, supportright(s)])
end
