# Implementation of
#  "Table-driven Implementation of the Logarithm Function in IEEE Floating-point Arithmetic"
#  Tang, Ping-Tak Peter
#  ACM Trans. Math. Softw. (1990), 16(4):378--400
#  http://dx.doi.org/10.1145/98267.98294

# Other implementations:
#  Apple OS X libm (APSL v2)
#   http://www.opensource.apple.com/source/Libm/Libm-2026/Source/Intel/xmm_log.c
#  libacml_mv (GPL v2)
#   http://svn.open64.net/svnroot/open64/trunk/osprey/libacml_mv/src/log1p.c
#  libclc (BSD): modified to use different table size
#   https://github.com/llvm-mirror/libclc/blob/master/generic/lib/math/log1p.cl

const t_log64 = Array((Float64,Float64),129)
N=39 # (can be up to N=42, which appears to be what Apple's libm uses).
sN = 2.0^N
isN = 1.0/sN
s7 = 2.0^7
is7 = 1.0/s7

for j=0:128
    l_big = Base.log(big(1.0+j*is7))
    l_hi = isN*float64(round(sN*l_big))
    l_lo = float64(l_big-l_hi)
    t_log64[j+1] = (l_hi,l_lo)
    # c_invF[j+1] = 1.0+is7*j
end

const t_log32 = Array((Float32,Float32),129)
N=16
sN = 2f0^N
isN = 1f0/sN
for j=0:128
    l_big = Base.log(big(1.0+j*is7))
    l_hi = isN*float32(round(sN*l_big))
    l_lo = float32(l_big-l_hi)
    t_log32[j+1] = (l_hi,l_lo)
end

# Procedure 2
@inline function log_proc2(f::Float64)
    ## Step 1
    g = 1.0/(2.0+f)
    u = 2.0*f*g
    v = u*u

    ## Step 2
    q = u*v*@horner(v,
                    0x1.5_5555_5555_54e6p-4,
                    0x1.9_9999_99ba_c6d4p-7,
                    0x1.2_4923_07f1_519fp-9,
                    0x1.c_8034_c85d_fff0p-12)

    ## Step 3
    # could improve this with an fma
    # return u+muladd(fma(-u,f,2(f-u)),g,q)
    u1 = float64(float32(u)) # round to 24 bits
    f1 = float64(float32(f))
    f2 = f-f1
    u2 = ((2.0*(f-u1)-u1*f1)-u1*f2)*g

    ## Step 4
    u1 + (u2 + q)
end

@inline function log_proc2(f::Float32)
    ## Step 1
    # compute in higher precision
    u64 = Float64(2f0*f)/(2.0+f)
    u = Float32(u64)
    v = u*u

    ## Step 2
    q = u*v*@horner(v,
                    Float32(0x1.555552p-4),
                    Float32(0x1.9a012ap-7))

    ## Step 3
    # unnecessary, as we have Float64 value

    ## Step 4
    Float32(u64 + q)
end

# Procedure 1
@inline function log_proc1(y::Float64,mf::Float64,F::Float64,f::Float64,jp::Int)
    ## Step 2
    @inbounds hi,lo = libm.t_log64[jp]
    l_hi = mf* 0x1.62e42fefa4p-1 + hi
    l_lo = mf*-0x1.8432a1b0e2634p-43 + lo

    ## Step 3
    # @inbounds u = f*c_invF[jp]
    # u = f/F
    # q = u*u*@horner(u,
    #                 -0x1.0_0000_0000_0001p-1,
    #                 +0x1.5_5555_5550_9ba5p-2,
    #                 -0x1.f_ffff_ffeb_6526p-3,
    #                 +0x1.9_99b4_dfed_6fe4p-3,
    #                 -0x1.5_5576_6647_2e04p-3)

    ## Step 3' (alternative)
    u = (2.0f)/(y+F)
    v = u*u
    q = u*v*@horner(v,
                    0x1.5_5555_5555_0286p-4,
                    0x1.9_99a0_bc71_2416p-7)

    ## Step 4
    l_hi + (u + (q + l_lo))
end

@inline function log_proc1(y::Float32,mf::Float32,F::Float32,f::Float32,jp::Int)
    # Procedure 1
    ## Step 2
    @inbounds hi,lo = t_log32[jp]
    l_hi = mf*Float32(0x1.62e4p-1) + hi
    l_lo = mf*Float32(0x1.7f7d1cp-20) + lo

    ## Step 3
    # @inbounds u = f*c_invF[jp]
    # q = u*u*@horner(u,
    #                 Float32(-0x1.00006p-1),
    #                 Float32(0x1.55546cp-2))

    ## Step 3' (alternative)
    u = (2f0f)/(y+F)
    v = u*u
    q = u*v*Float32(0x1.555584p-4)

    ## Step 4
    l_hi + (u + (q + l_lo))
end


function log_tang(x::Float64)
    if x > 0.0
        x == Inf && return x

        # Step 2
        if 0x1.e_0fab_fbc7_02a3p-1 < x < 0x1.1_082b_577d_34eep0
            f = x-1.0
            return log_proc2(f)
        end

        # Step 3
        xu = reinterpret(Uint64,x)
        m = int(xu >> 52) & 0x07ff
        if m == 0 # x is subnormal
            x *= 0x1p54 # normalise significand
            xu = reinterpret(Uint64,x)
            m = int(xu >> 52) & 0x07ff - 54
        end
        m -= 1023
        y = reinterpret(Float64,(xu & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000)

        mf = Float64(m)
        F = (y + 0x1p45) - 0x1p45 # 0x1p-7*round(0x1p7*y)
        f = y-F
        jp = itrunc(0x1p7*F)-127

        return log_proc1(y,mf,F,f,jp)
    elseif x == 0.0
        -Inf
    elseif isnan(x)
        NaN
    else
        throw(DomainError())
    end
end

function log_tang(x::Float32)
    if x > 0f0
        x == Inf32 && return x

        # Step 2
        if Float32(0x1.e0fabep-1) < x < Float32(0x1.1082b6p+0)
            f = x-1f0
            return log_proc2(f)
        end

        # Step 3
        xu = reinterpret(Uint32,x)
        m = int(xu >> 23) & 0x00ff
        if m == 0 # x is subnormal
            x *= Float32(0x1p25) # normalise significand
            xu = reinterpret(Uint32,x)
            m = int(xu >> 23) & 0x00ff - 25
        end
        m -= 127
        y = reinterpret(Float32,(xu & 0x007f_ffff) | 0x3f80_0000)

        mf = Float32(m)
        F = (y + Float32(0x1p16)) - Float32(0x1p16) # 0x1p-7*round(0x1p7*y)
        f = y-F
        jp = itrunc(Float32(0x1p7)*F)-127

        log_proc1(y,mf,F,f,jp)
    elseif x == 0f0
        -Inf32
    elseif isnan(x)
        NaN32
    else
        throw(DomainError())
    end
end

function log1p_tang(x::Float64)
    if x > -1.0
        x == Inf && return x
        if -0x1p-53 < x < 0x1p-53
            return x # Inexact

        # Step 2
        elseif -0x1.f_0540_438f_d5c4p-5 < x < 0x1.0_82b5_77d3_4ed8p-4
            return log_proc2(x)
        end

        # Step 3
        z = 1.0 + x
        zu = reinterpret(Uint64,z)
        m = int(zu >> 52) & 0x07ff - 1023 # z cannot be subnormal
        c = m > 0 ? 1.0-(z-x) : x-(z-1.0) # 1+x = z+c
        y = reinterpret(Float64,(zu & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000)

        mf = Float64(m)
        F = (y + 0x1p45) - 0x1p45 # 0x1p-7*round(0x1p7*y)
        f = (y - F) + ldexp(c,-m) #2^m(F+f) = 1+x = z+c
        jp = itrunc(0x1p7*F)-127

        log_proc1(y,mf,F,f,jp)
    elseif x == -1.0
        -Inf
    elseif isnan(x)
        NaN
    else
        throw(DomainError())
    end
end

function log1p_tang(x::Float32)
    if x > -1f0
        x == Inf32 && return x
        if -0x1p-53 < x < 0x1p-53
            return x # Inexact

        # Step 2
        elseif Float32(-0x1.f05406p-5) < x < Float32(0x1.082b58p-4)
            return log_proc2(x)
        end

        # Step 3
        z = 1f0 + x
        zu = reinterpret(Uint32,z)
        m = int(zu >> 23) & 0x00ff - 127 # z cannot be subnormal
        c = m > 0 ? 1f0-(z-x) : x-(z-1f0) # 1+x = z+c
        y = reinterpret(Float32,(xu & 0x007f_ffff) | 0x3f80_0000)

        mf = Float32(m)
        F = (y + Float32(0x1p16)) - Float32(0x1p16) # 0x1p-7*round(0x1p7*y)
        f = (y - F) + ldexp(c,-m) #2^m(F+f) = 1+x = z+c
        jp = itrunc(Float32(0x1p7)*F)-127

        log_proc1(y,mf,F,f,jp)
    elseif x == -1f0
        -Inf32
    elseif isnan(x)
        NaN32
    else
        throw(DomainError())
    end
end

@vectorize_1arg Real log_tang
@vectorize_1arg Real log1p_tang
