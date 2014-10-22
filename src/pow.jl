# a rough translation of openlibm/src/e_pow.c under following licence:

## Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
##
## Developed at SunSoft, a Sun Microsystems, Inc. business.
## Permission to use, copy, modify, and distribute this
## software is freely granted, provided that this notice 
## is preserved.

# note: this doesn't handle edge cases or negative numbers correctly.

function log2_ext(x)
    sig,n = frexp(x)
    ax = 2sig; n -= 1

    # reduce to 
    # k = 0, sqrt(3)/2 = 0.866 < ax < sqrt(3/2) = 1.22
    # k = 1, sqrt(3/2) = 1.22 < ax < sqrt(3) = 1.73

    bp = 1.0
    dp_h = 0.0
    dp_l = 0.0

    if ax < sqrt(3/2)
    elseif ax < sqrt(3)
        bp = 1.5
        # log2(1.5)
        dp_h = 5.84962487220764160156e-01
        dp_l = 1.35003920212974897128e-08
    else
        n += 1
        ax *= 0.5
    end
    # log(x) = 2s + 2/3 s^3 + 2/5 s^5...

    u = ax - bp # usual f, reduced to -0.28 < u < 0.23
    v = 1.0/(ax + bp) 
    ss = u*v # ss = (x-1)/(x+1) or (x-1.5)/(x+1.5) < 0.125 = 2^-3 
    s_h = trunc32(ss)
    t_h = trunc32(ax+bp) # according to comment
    t_l = ax - (t_h-bp)
    s_l = v*((u-s_h*t_h)-s_h*t_l) # division correction?

    z = ss*ss # < 2^-6
    # log(ax) = 2s + 2/3*(s^3 + r) = 2/3s (3 + s^2 + r) 
    r = z*z*@horner(z,
                    5.99999999999994648725e-01,
                    4.28571428578550184252e-01,
                    3.33333329818377432918e-01,
                    2.72728123808534006489e-01,
                    2.30660745775561754067e-01,
                    2.06975017800338417784e-01) # < 2^-12
    r += s_l*(s_h+ss)
    s2 = s_h*s_h
    # t = @split_add(3.0,s2,r)
    t_h = trunc32(3.0 + s2 + r)
    t_l = r+(s2+(3.0-t_h))
    # u+v = ss*(1+...)
    # p = @split t*s
    u = s_h*t_h
    v = s_l*t_h+t_l*ss

    # 2/(3log2)
    cp    =  9.61796693925975554329e-01 # /* 0x3FEEC709, 0xDC3A03FD =2/(3ln2) */
    cp_h  =  9.61796700954437255859e-01 # /* 0x3FEEC709, 0xE0000000 =(float)cp */
    cp_l  = -7.02846165095275826516e-09 # /* 0xBE3E2FE0, 0x145B01F5 =tail of cp_h*/
        
    # 2/(3log2)*(ss+...)
    p_h = trunc32(u+v)
    p_l = v-(p_h-u)
    z_h = cp_h*p_h 
    z_l = cp_l*p_h + p_l*cp + dp_l

    t = float64(n)
    t1 = trunc32(((z_h+z_l) + dp_h) + t)
    t2 = z_l - (((t1 - t) - dp_h) - z_h)
    Double64(t1,t2)
end


function mul_ext(t::Double64,y::Float64)
    # split y, y*(t1+t2)
    y1 = trunc32(y)
    p_l = (y-y1)*t.hi + y*t.lo
    p_h = y1*t.hi
    Double64(p_h,p_l)
end

function exp2(p::Double64)
    z = p.lo + p.hi
    # check for under/overflow

    # 2^(p_l+p_h)

    zr = round(z)
    n = int(zr)
    p_h = p.hi - zr
    p_l = p.lo

    lg2  =  6.93147180559945286227e-01 # /* 0x3FE62E42, 0xFEFA39EF */
    lg2_h  =  6.93147182464599609375e-01 # /* 0x3FE62E43, 0x00000000 */
    lg2_l  = -1.90465429995776804525e-09 # /* 0xBE205C61, 0x0CA86C39 */


    # 1) t1 = 2(e^z-1-z)/(e^z-1) = z - z^2/6 + z^4/360 - z^6/15120 + ...
    # 2) r = z*t1/(t1-2) = 1 - e^z + z = -(z^2/2 + s^3/6 + ...)
    # 3) e^z = 1 - r + z
    t = trunc32(p_l + p_h)
    u = t*lg2_h
    v = (p_l-(t-p_h))*lg2+t*lg2_l
    z = u+v
    w = v-(z-u);
    t  = z*z;
    t1  = z - t*@horner(t,
                        1.66666666666666019037e-01,
                        -2.77777777770155933842e-03,
                        6.61375632143793436117e-05,
                        -1.65339022054652515390e-06,
                        4.13813679705723846039e-08)
    r = (z*t1)/(t1 - 2.0) - (w+z*w)
    z = 1.0-(r-z)
    scalbn(z,n)
end

pow(x,y) = exp2(mul_ext(log2_ext(x),y))

