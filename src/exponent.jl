
function exponent(x::Float64)
    xu = reinterpret(Uint64,x)
    int(xu >> 52) & 0x07ff - 1023
end

function exponent(x::Float32)
    xu = reinterpret(Uint32,x)
    int(xu >> 23) & 0x00ff -126
end

function setexponent(x::Float64,k::Int)
    xu = reinterpret(Uint64,x)
    ku = uint64((k + 1023) & 0x07ff) << 52
    xu = (xu & 0x800f_ffff_ffff_ffff) | ku
    reinterpret(Float64,xu)
end
