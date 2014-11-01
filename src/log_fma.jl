function log_fma(x::Float64)
    if x > 0.0
        if isinf(x) 
            Inf
        else
            xu = reinterpret(Uint64,x)
            k = int(xu >> 52) & 0x07ff
            if k == 0 # x is subnormal
                x *= 0x1p54 # normalise significand
                xu = reinterpret(Uint64,x)
                k = int(xu >> 52) & 0x07ff - 54
            end
            k -= 1022
            xu = (xu & 0x800f_ffff_ffff_ffff) | 0x3fe0_0000_0000_0000
            s = reinterpret(Float64,xu)
            if s < 0.7071067811865476
                s *= 2.0
                k -= 1
            end
            f = s-1.0
            fk = float64(k)

            ln2_hi = 6.93147180369123816490e-01
            ln2_lo = 1.90821492927058770002e-10

            lg1 = 6.666666666666735130e-01
            lg2 = 3.999999999940941908e-01
            lg3 = 2.857142874366239149e-01
            lg4 = 2.222219843214978396e-01
            lg5 = 1.818357216161805012e-01
            lg6 = 1.531383769920937332e-01
            lg7 = 1.479819860511658591e-01

            hfsq = 0.5*f*f
            s = f/(2.0+f)
            z = s*s
            w = z*z
            t1 = w*@horner_fma(w,lg2,lg4,lg6)
            t2 = @horner_fma(w,lg1,lg3,lg5,lg7)
            R = fma(t2,z,t1)
            fk*ln2_hi+(f-(hfsq-(s*(hfsq+R)+(fk*ln2_lo))))
        end
    elseif x == 0.0
        -Inf
    else
        throw(DomainError())
    end
end

function log_muladd(x::Float64)
    if x > 0.0
        if isinf(x) 
            Inf
        else
            xu = reinterpret(Uint64,x)
            k = int(xu >> 52) & 0x07ff
            if k == 0 # x is subnormal
                x *= 0x1p54 # normalise significand
                xu = reinterpret(Uint64,x)
                k = int(xu >> 52) & 0x07ff - 54
            end
            k -= 1022
            xu = (xu & 0x800f_ffff_ffff_ffff) | 0x3fe0_0000_0000_0000
            s = reinterpret(Float64,xu)
            if s < 0.7071067811865476
                s *= 2.0
                k -= 1
            end
            f = s-1.0
            fk = float64(k)

            ln2_hi = 6.93147180369123816490e-01
            ln2_lo = 1.90821492927058770002e-10

            lg1 = 6.666666666666735130e-01
            lg2 = 3.999999999940941908e-01
            lg3 = 2.857142874366239149e-01
            lg4 = 2.222219843214978396e-01
            lg5 = 1.818357216161805012e-01
            lg6 = 1.531383769920937332e-01
            lg7 = 1.479819860511658591e-01

            hfsq = 0.5*f*f
            s = f/(2.0+f)
            z = s*s
            w = z*z
            t1 = w*@horner_muladd(w,lg2,lg4,lg6)
            t2 = @horner_muladd(w,lg1,lg3,lg5,lg7)
            R = muladd(t2,z,t1)
            fk*ln2_hi+(f-(hfsq-(s*(hfsq+R)+(fk*ln2_lo))))
        end
    elseif x == 0.0
        -Inf
    else
        throw(DomainError())
    end
end
