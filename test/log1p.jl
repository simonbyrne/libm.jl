using libm

N = 1_000_000
srand(1)
X = exp(randn(N))
function maxerr{T}(X::AbstractArray{T})
    e = 0.0
    for x in X
        a = libm.log1p_tang(x)
        b = log1p(big(x))
        d = abs(T(b-a))
        r = d / eps(T(b))
        e = max(e,r)
    end
    e
end

@printf "Float64 max rel. error %f ulps\n" maxerr(X)
@printf "Float32 max rel. error %f ulps\n" maxerr(float32(X))
