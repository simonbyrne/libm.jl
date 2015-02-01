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
        @time test($X)
        println()
    end
end

syslog1p(x::Float64) = Core.Intrinsics.nan_dom_err(ccall((:log1p,:libm),Float64,(Float64,),x), x)
syslog1p(x::Float32) = Core.Intrinsics.nan_dom_err(ccall((:log1pf,:libm),Float32,(Float32,),x), x)

println("Float64")
println("-------------------")
@testsum "Openlibm:" log1p X
@testsum "System:" syslog1p X

#@testsum "Julia fdlibm:" libm.log X
@testsum "Julia tang:" libm.log1p_tang X
# @testsum "muladd:" libm.log_muladd X

println()
println()

X = float32(X)

println("Float32")
println("-------------------")
@testsum "Openlibm:" log1p X
@testsum "System:" syslog1p X

#@testsum "Julia fdlibm:" libm.log X
@testsum "Julia tang:" libm.log1p_tang X
