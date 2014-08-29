
tiny(x::Float64) = 0x1p-1023
huge(x::Float64) = 0x1p1023


# scalbn(x,n) = x * 2^n
function scalbn(x::Float64, n::Int)
    xu = reinterpret(Uint64,x)
    k = int(xu >> 52) & 0x07ff
    if k == 0 # x is subnormal
        x == zero(x) && return x
        x *= 0x1p54 # normalise significand
        xu = reinterpret(Uint64,x)
        k = int(xu >> 52) & 0x07ff
        k -= 54
        n < -50000 && return tiny(x)*x # so exponent doesn't underflow
    end
    k == 0x07ff && return x + x # NaN or Inf
    k += n    
    if k > 0x07fe || n > 50000  # overflow result
        return huge(x)*copysign(huge(x),x)
    elseif k > 0 # normal result
        xu = (xu & 0x800f_ffff_ffff_ffff) | (int64(k) << 52)
        return reinterpret(Float64,xu)
    elseif k > -54 # subnormal result
        xu = (xu & 0x800f_ffff_ffff_ffff) | (int64(k+54) << 52)
        return reinterpret(Float64,xu)*0x1p-54
    else # underlow result
        return tiny(x)*copysign(tiny(x),x)
    end
end
