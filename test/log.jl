using libm

N = 1_000_000
srand(1)
X = exp(randn(N))
function maxerr(X)
    e = 0.0
    for x in X
        a = libm.log_tang(x)
        b = log(big(x))
        d = abs(Float64(b-a))
        r = d / eps(Float64(b))
        e = max(e,r)
    end
    e
end

@printf "Max rel. error %f ulps\n" maxerr(X)
