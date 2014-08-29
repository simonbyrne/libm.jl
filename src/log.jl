# a rough translation of openlibm/src/e_log.c under following licence:

## Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
##
## Developed at SunSoft, a Sun Microsystems, Inc. business.
## Permission to use, copy, modify, and distribute this
## software is freely granted, provided that this notice 
## is preserved.


# reduce to x = 2^k(1+f) where sqrt(2)/2 < 1+f < sqrt(2)
function log_red(x::Float64)
    s, k = frexp(x)
    if s < 0.7071067811865476
        s *= 2.0
        k -= 1
    end
    k, s-1.0
end

# compute log of reduced quantity
function log1p_kernel(k::Int,f::Float64)
    ln2_hi = 6.93147180369123816490e-01
    ln2_lo = 1.90821492927058770002e-10

    hfsq = 0.5*f*f
    s = f/(2.0+f)
    z = s*s
    R = z*@horner(z,
                6.666666666666735130e-01,
                3.999999999940941908e-01,
                2.857142874366239149e-01,
                2.222219843214978396e-01,
                1.818357216161805012e-01,
                1.531383769920937332e-01,
                1.479819860511658591e-01)
    k*ln2_hi+(f-(hfsq-(s*(hfsq+R)+(k*ln2_lo))))
end

# function
function log(x::Float64)
    if x > 0.0
        if isinf(x) 
            Inf
        else
            k,f = log_red(x)
            log1p_kernel(k,f)
        end
    elseif x == 0.0
        -Inf
    else
        NaN
    end
end
        
