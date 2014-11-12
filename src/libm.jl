module libm

import Base.Math.@horner

macro horner_fma(x, p...)
    ex = esc(p[end])
    for i = length(p)-1:-1:1
        ex = :(fma(t,$ex,$(esc(p[i]))))
    end
    Expr(:block, :(t = $(esc(x))), ex)
end

macro horner_muladd(x, p...)
    ex = esc(p[end])
    for i = length(p)-1:-1:1
        ex = :(muladd(t,$ex,$(esc(p[i]))))
    end
    Expr(:block, :(t = $(esc(x))), ex)
end


half(x::Float64) = 0.5
half(x::Float32) = 0.5f0

immutable Double{T}
    hi::Float64
    lo::Float64
end

typealias Double64 Double{Float64}

trunc32(x::Float64) = reinterpret(Float64,reinterpret(Uint64,x) & 0xffff_ffff_0000_0000)

# include("frexp.jl")
# include("exponent.jl")
# include("scalbn.jl")

# include("trunc.jl")
# include("round.jl")

include("log.jl")
# include("log_fma.jl")
include("log_tang.jl")
include("expm1.jl")

include("pow.jl")
include("pow1p.jl")

end # module
