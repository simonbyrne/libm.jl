function frexp(x::Float64)
    xu = reinterpret(Uint64,x)
    k = int(xu >> 52) & 0x07ff
    if k == 0 # x is subnormal
        x == zero(x) && return x,0
        x *= 0x1p54 # normalise significand
        xu = reinterpret(Uint64,x)
        k = int(xu >> 52) & 0x07ff - 54
    elseif k == 0x07ff
        # NaN or Inf
        return x,0
    end
    k -= 1022
    xu = (xu & 0x800f_ffff_ffff_ffff) | 0x3fe0_0000_0000_0000
    reinterpret(Float64,xu), k
end

function frexp(x::Float32)
    xu = reinterpret(Uint32,x)
    k = int(xu >> 23) & 0x00ff
    if k == 0 # x is subnormal
        x == zero(x) && return x,0
        x *= float32("0x1p25") # no Float32 hex literal
        xu = reinterpret(Uint32,x)
        k = int(xu >> 23) & 0x00ff - 25
    elseif k == 0x00ff
        # NaN or Inf
        return x,0
    end
    k -= 126
    xu = (xu & 0x807f_ffff) | 0x3f00_0000
    reinterpret(Float32,xu), k
end


