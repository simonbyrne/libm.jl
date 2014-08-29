module libm

import Base.Math.@horner

half(x::Float64) = 0.5
half(x::Float32) = 0.5f0


include("frexp.jl")
include("exponent.jl")
include("scalbn.jl")

include("trunc.jl")
include("round.jl")

include("log.jl")

include("pow.jl")

end # module
