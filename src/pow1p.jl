
function log21p_ext(x)
    u = 1.0 + x
    sig,n = frexp(u)
    au = 2sig; n -= 1

    # reduce to 
    # k = 0, sqrt(3)/2 = 0.866 < au < sqrt(3/2) = 1.22
    # k = 1, sqrt(3/2) = 1.22 < au < sqrt(3) = 1.73

    if au < sqrt(3/2)
        bp = 1.0
        dp_h = 0.0
        dp_l = 0.0
    elseif au < sqrt(3)
        bp = 1.5
        # log2(1.5)
        dp_h = 5.84962487220764160156e-01
        dp_l = 1.35003920212974897128e-08
    else
        bp = 1.0
        dp_h = 0.0
        dp_l = 0.0

        n += 1
        au *= 0.5
    end

    # exact 1.0 + x = u + c
    if x >= 1.0 
        # u <= 2x
        c = 1.0-(u-x)
    else
        # 1) x > -0.5: 0.5 <= u <= 2.0
        # 2) x <= -0.5: u = x+1.0 exact, hence so is u-1.0
        c = x-(u-1.0) 
    end

    # x < -0.5: c = 0
    # x < 0.5: |c| < eps(1.0)
    # x < 0x1p53: |c| in {-eps(x), 0, eps(x)} 
    # otherwise: c = 1

    # scale to au:  1 + x = (au + ac)*2^n
    ac = ldexp(c,-n) 
    # log(x) = 2s + 2/3 s^3 + 2/5 s^5...
    f = au - bp # usual f, reduced to -0.28 < u < 0.23
    f1 = f + ac 
    f2 = ac + (f-f1)
    g = au + bp
    g_h = trunc32(g) # according to comment
    g_l = ac + (au + (bp - g_h))
    rg = 1.0/g

    ss = f1*rg # ss = (x-1)/(x+1) or (x-1.5)/(x+1.5) < 0.125 = 2^-3 
    s_h = trunc32(ss)
    s_l = (((f1-s_h*g_h)-s_h*g_l) + f2)*rg # division correction?

    z = ss*ss # < 2^-6
    # log(au) = 2s + 2/3*(s^3 + r) = 2/3s (3 + s^2 + r) 
    r = z*z*@horner(z,
                    5.99999999999994648725e-01,
                    4.28571428578550184252e-01,
                    3.33333329818377432918e-01,
                    2.72728123808534006489e-01,
                    2.30660745775561754067e-01,
                    2.06975017800338417784e-01) # < 2^-12
    r += s_l*(s_h+ss)
    s2 = s_h*s_h
    t_h = 3.0 + s2 + r
    t_h = trunc32(t_h)
    t_l = r-((t_h-3.0)-s2)
    # u+v = ss*(1+...)
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

function exp2m1(p::Double64)
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
    (scalbn(1.0,n) - 1.0) - scalbn(r-z,n)
end


pow1p(x,y) = exp2(mul_ext(log21p_ext(x),y))
powm1(x,y) = exp2m1(mul_ext(log2_ext(x),y))
pow1pm1(x,y) = exp2m1(mul_ext(log21p_ext(x),y))
