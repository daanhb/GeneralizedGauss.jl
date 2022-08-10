
function oscillatory_basis_exp(omega, n)
    x = IdentityFunction()
    w1 = exp(im*omega*x)
    w2 = exp(-im*omega*x)
    if iseven(n)
        w1 * ChebyshevT(n>>1) ⊕ w2 * ChebyshevT(n>>1)
    else
        w1 * ChebyshevT(n>>1+1) ⊕ w2 * ChebyshevT(n>>1)
    end
end

function oscillatory_basis_cos(omega, n)
    x = IdentityFunction()
    w1 = cos(omega*x)
    w2 = sin(omega*x)
    if iseven(n)
        w1 * ChebyshevT(n>>1) ⊕ w2 * ChebyshevT(n>>1)
    else
        w1 * ChebyshevT(n>>1+1) ⊕ w2 * ChebyshevT(n>>1)
    end
end

function log_basis(n)
    b1 = ChebyshevT(n>>1) → 0..1
    b2 = ChebyshevT(n>>1+1) → 0..1
    iseven(n) ? b1 ⊕ Log() * b1 : b2 ⊕ Log() * b1
end

function exp_basis(n)
    b1 = Chebyhev(n>>1) → 0..1
    b2 = Chebyhev(n>>1+1) → 0..1
    iseven(n) ? b1 ⊕ Exp() * b1 : b2 ⊕ Exp() * b1
end
