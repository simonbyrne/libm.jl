# A modification of Tang's table-driven logarithm for different base.

const t_log2_64 = Array((Float64,Float64),129)
N=39 # (can be up to N=42, which appears to be what Apple's libm uses).
sN = 2.0^N
isN = 1.0/sN
s7 = 2.0^7
is7 = 1.0/s7

for j=0:128
    l_big = Base.log2(big(1.0+j*is7))
    l_hi = isN*float64(round(sN*l_big))
    l_lo = float64(l_big-l_hi)
    t_log2_64[j+1] = (l_hi,l_lo)
    # c_invF[j+1] = 1.0+is7*j
end

const inv_ln2 = 0x1.71547652b82fep+0
const inv_ln2_hi = 0x1.7154766p+0 # 28 bits
const inv_ln2_lo = -0x1.a8fa03d1105ep-29

# const t_log_32 = Array((Float32,Float32),129)
# N=16
# sN = 2f0^N
# isN = 1f0/sN
# for j=0:128
#     l_big = Base.log(big(1.0+j*is7))
#     l_hi = isN*float32(round(sN*l_big))
#     l_lo = float32(l_big-l_hi)
#     t_log_32[j+1] = (l_hi,l_lo)
# end

# Procedure 1
@inline function log2_proc1(y::Float64,mf::Float64,F::Float64,f::Float64,jp::Int)
    ## Steps 1 and 2
    @inbounds hi,lo = libm.t_log2_64[jp]
    l_hi = mf + hi
    l_lo = lo

    ## Step 3' (alternative)
    u = (2.0f)/(y+F)
    us = inv_ln2*u
    v = u*u
    q = us*v*@horner(v,
                    0x1.5_5555_5555_0286p-4,
                    0x1.9_99a0_bc71_2416p-7)

    ## Step 4
    l_hi + (us + (q + l_lo))
end

@inline function log2_proc2(f::Float64)
    ## Step 1
    g = 1.0/(2.0+f)
    u = 2.0*f*g
    us = 1.4426950408889634*u
    v = u*u

    ## Step 2
    q = us*v*@horner(v,
                    0x1.5_5555_5555_54e6p-4,
                    0x1.9_9999_99ba_c6d4p-7,
                    0x1.2_4923_07f1_519fp-9,
                    0x1.c_8034_c85d_fff0p-12)

    ## Step 3
    u1 = float64(float32(u)) # round to 24 bits
    f1 = float64(float32(f))
    f2 = f-f1
    u2 = ((2.0*(f-u1)-u1*f1)-u1*f2)*g

    ## Step 4
    u1*inv_ln2_hi + (u2*inv_ln2 + (u1*inv_ln2_lo + q))
end


function log2_tang(x::Float64)
    if x > 0.0
        x == Inf && return x

        # Step 2
        if 0x1.e_0fab_fbc7_02a3p-1 < x < 0x1.1_082b_577d_34eep0
            f = x-1.0
            return log2_proc2(f)
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

        log2_proc1(y,mf,F,f,jp)
    elseif x == 0.0
        -Inf
    elseif isnan(x)
        NaN
    else
        throw(DomainError())
    end
end

@vectorize_1arg Real log2_tang
