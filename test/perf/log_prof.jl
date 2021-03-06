using libm
srand(1)

X = exp(randn(10_000_000))

macro testsum(name,fn,X)
    quote
        function test(X)
            s = 0.0
            for x in X
                s += $fn(x)
            end
            s
        end

        test($X)
        test($X)

        println($name)
        @profile test($X)
        Profile.print()
    end
end

# syslog(x::Float64) = Core.Intrinsics.nan_dom_err(ccall((:log,:libm),Float64,(Float64,),x), x)

# @testsum "Openlibm:" log X
# @testsum "System:" syslog X
# @testsum "Julia fdlibm:" libm.log X
@testsum "Julia tang:" libm.log_tang X
# @testsum "fma:" libm.log_fma X
# @testsum "muladd:" libm.log_muladd X
