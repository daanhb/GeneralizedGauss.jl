
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
    if iseven(n)
        ChebyshevT(n>>1, 0, 1) ⊕ Log() * ChebyshevT(n>>1, 0, 1)
    else
        ChebyshevT(n>>1+1, 0, 1) ⊕ Log() * ChebyshevT(n>>1, 0, 1)
    end
end

function exp_basis(n)
    if iseven(n)
        ChebyshevT(n>>1, 0, 1) ⊕ Exp() * ChebyshevT(n>>1, 0, 1)
    else
        ChebyshevT(n>>1+1, 0, 1) ⊕ Exp() * ChebyshevT(n>>1, 0, 1)
    end
end
