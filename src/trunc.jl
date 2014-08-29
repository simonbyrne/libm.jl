function trunc(x::Float64)
    xu = reinterpret(Uint64,x)
    k = int(xu >> 52) & 0x07ff - 1023
    if k > 1023 # Inf or NaN
        x+x
    elseif k < 0
        copysign(0.0,x)
    else
        reinterpret(Float64, xu & ~(0x000f_ffff_ffff_ffff >> k))
    end
end
        
function trunc(x::Float32)
    xu = reinterpret(Uint32,x)
    k = int(xu >> 23) & 0x00ff - 127
    if k > 127 # Inf or NaN
        x+x
    elseif k < 0
        copysign(0f0,x)
    else
        reinterpret(Float32, xu & ~(0x007f_ffff >> k))
    end
end
        
