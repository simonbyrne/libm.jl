function round(x::FloatingPoint)
    isfinite(x) || return x+x
    t = trunc(x)
    if x > 0.0
        x-t >= half(x) ? t+one(t) : t
    else
        x-t <= -half(x) ? t-one(t) : t
    end
end

function floor(x::FloatingPoint)
    t = trunc(x)
    x < zero(x) && x!=t ? t-one(x) : t
end


function ceil(x::FloatingPoint)
    t = trunc(x)
    x > zero(x) && x!=t ? t+one(x) : t
end
