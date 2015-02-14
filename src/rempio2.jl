
immutable Double64
    hi::Float64
    lo::Float64
end

function rem_pio2(x::Float64)
    two24 =  1.67772160000000000000e+07   # 0x4170_0000_0000_0000
    invpio2 =  6.36619772367581382433e-01 # 0x3FE4_5F30_6DC9_C883
    pio2_1  =  1.57079632673412561417e+00 # 0x3FF9_21FB_5440_0000
    pio2_1t =  6.07710050650619224932e-11 # 0x3DD0_B461_1A62_6331
    pio2_2  =  6.07710050630396597660e-11 # 0x3DD0_B461_1A60_0000
    pio2_2t =  2.02226624879595063154e-21 # 0x3BA3_198A_2E03_7073
    pio2_3  =  2.02226624871116645580e-21 # 0x3BA3_198A_2E00_0000
    pio2_3t =  8.47842766036889956997e-32 # 0x397B_839A_2520_49C1

    ux = reinterpret(UInt64,x)
    hx = (ux >> 32) % UInt
    ix = hx & 0x7fff_ffff
    
    if ix <= 0x3fe9_21fb # |x| ~<= pi/4
        return 0, Double64(x, 0.0)
    elseif ix <= 0x400f_6a7a # |x| ~<= 5pi/4
        if ((ix & 0xf_ffff) == 0x9_21fb) # |x| ~= pi/2, pi/4 
            @goto medium
        elseif ix <= 0x4002_d97c # |x| ~<= 3pi/4
            if hx > 0
                z = x - pio2_1
                y_hi = z - pio2_1t
                y_lo = (z - y_hi) - pio2_1t
                return 1, Double64(y_hi, y_lo)
            else
                z = x + pio2_1
                y_hi = z + pio2_1t
                y_lo = (z - y_hi) + pio2_1t
                return -1, Double64(y_hi, y_lo)
            end
        else # 3pi/4 ~<= |x| ~<= 5pi/4
            if hx > 0
                z = x - 2pio2_1
                y_hi = z - 2pio2_1t
                y_lo = (z - y_hi) - 2pio2_1t
                return 2, Double64(y_hi,y_lo)
            else
                z = x + 2pio2_1
                y_hi = z + 2pio2_1t
                y_lo = (z - y_hi) + 2pio2_1t
                return -2, Double64(y_hi,y_lo)
            end
        end
    elseif ix <= 0x401c_463b # |x| ~<= 9pi/4
        if ix <= 0x4015_fdbc # |x| ~<= 7pi/4
            if ix == 0x4012_d97c # |x| ~= 3pi/2
                @goto medium
            elseif hx > 0
                z = x - 3pio2_1
                y_hi = z - 3pio2_1t
                y_lo = (z - y_hi) - 3pio2_1t
                return 3, Double64(y_hi,y_lo)
            else
                z = x + 3pio2_1
                y_hi = z + 3pio2_1t
                y_lo = (z - y_hi) + 3pio2_1t
                return -3, Double64(y_hi,y_lo)
            end
        else # 7pi/4 ~<= |x| ~<= 9pi/4
            if ix == 0x4019_21fb # |x| ~= 4pi/2
                @goto medium
            elseif hx > 0
                z = x - 4pio2_1
                y_hi = z - 4pio2_1t
                y_lo = (z - y_hi) - 4pio2_1t
                return 4, Double64(y_hi,y_lo)
            else
                z = x + 4pio2_1
                y_hi = z + 4pio2_1t
                y_lo = (z - y_hi) + 4pio2_1t
                return -4, Double64(y_hi,y_lo)
            end
        end
    end

    if ix < 0x413921fb # |x| ~< 2^20*(pi/2)
        @label medium

        fn = round(x*invpio2) # assumes fast, ties-to-even
        n = unsafe_trunc(Int,fn)
        r = x-fn*pio2_1
        w = fn*pio2_1t

        j = ix >> 20 # biased exponent of x
        y_hi = r-w
        k = ((reinterpret(UInt64,y_hi) >> 52) % UInt) & 0x7ff # biased exponent of y_hi
        i = j - k
        if i > 16 # require 2nd iteration
            t = r
            w = fn*pio2_2
            r = t-w
            w = fn*pio2_2t - ((t-r)-w)
            y_hi = r-w
            k = UInt(reinterpret(UInt64,y_hi) >> 52) & 0x7ff # biased exponent of y_hi
            i = j - k
            if i > 49 # require 3rd iteration
                t = r
                w = fn*pio2_3
                r = t-2
                w = fn*pio2_3t - ((t-r)-w)
                y_hi = r-w
            end
        end
        y_lo = (r-y_hi)-w
        return n, Double64(y_hi, y_lo)

    elseif ix >= 0x7ff00000 # Inf or NaN
        y = x-x
        return 0, Double64(y,y)

    else # TODO: slow method for large exponents
            
        return 0, Double64(0.0,0.0)
    end
end

