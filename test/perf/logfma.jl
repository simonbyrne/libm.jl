using libm

X = exp(randn(1_000_000))

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
        @time test($X)
        println()
    end
end

syslog(x::Float64) = ccall((:log,:libm),Float64,(Float64,),x)

@testsum "Openlibm:" log X
@testsum "System:" syslog X
@testsum "Julia:" libm.log X
@testsum "fma:" libm.log_fma X
@testsum "muladd:" libm.log_muladd X
