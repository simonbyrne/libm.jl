# "Table-driven Implementation of the Exponential Function in IEEE Floating-point Arithmetic"
# Tang, Ping-Tak Peter
# ACM Trans. Math. Softw., 1989, 15(2):144--157
# http://doi.acm.org/10.1145/63522.214389

# implemented in
#  Apple OS X libm (APSL v2)
#   http://www.opensource.apple.com/source/Libm/Libm-2026/Source/Intel/xmm_log.c
#  libacml_mv (GPL v2)
#   http://svn.open64.net/svnroot/open64/trunk/osprey/libacml_mv/src/log1p.c
#  libclc (BSD): modified to use different table size
#   https://github.com/llvm-mirror/libclc/blob/master/generic/lib/math/log1p.cl

const c_lead = Array(Float64,129)
const c_trail = Array(Float64,129)

N=39 # (can be up to N=42, which appears to be what Apple's libm uses).
sN = 2.0^N
s7 = 2.0^7
isN = 1.0/sN
is7 = 1.0/s7

for j=0:128
    l_big = Base.log(big(1.0+j*is7))
    l_hi = isN*float64(round(sN*l_big))
    c_lead[j+1] = l_hi
    c_trail[j+1] = float64(l_big-l_hi)
end

log2_big = Base.log(big(2.0))
const log2_lead = isN*float64(round(sN*log2_big))
const log2_trail = float64(log2_big - log2_lead)




function log_tang(x::Float64)
    if x > 0.0
        isinf(x) && return Inf

        # Step 2
        if 0x1.e_0fab_fbc7_02a3p-1 < x < 0x1.1_082b_577d_34eep0
            # Procedure 2
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
            # e.g. u2 = fma(-u1,f,2(f-u1))*g
            # or return fma(fma(-u,f,2(f-u)),g,q) ?
            u1 = float64(float32(u)) # round to 24 bits
            f1 = float64(float32(f))
            f2 = f-f1
            u2 = ((2.0*(f-u1)-u1*f1)-u1*f2)*g
            
            ## Step 4
            return u1 + (u2 + q)
        end
        
        # Step 3
        m = float(exponent(x))
        y = ldexp(x,-m) # y ∈ [1,2)
        F = (y + 0x1p45) - 0x1p45 # is7*round(s7*y)
        f = y-F
        jp = itrunc(s7*y)-127
        
        # Procedure 1
        ## Step 2
        l_lead = m*log2_lead + c_lead[jp]
        l_trail = m*log2_trail + c_trail[jp]

        ## Step 3
        # u = f/F # can tabulate 1/F
        # q = u*u*@horner(u,
        #                 +0x1.0_0000_0000_0001p-1,
        #                 +0x1.5_5555_5550_9ba5p-2,
        #                 -0x1.f_ffff_ffeb_6526p-3,
        #                 +0x1.9_99b4_dfed_6fe4p-3,
        #                 -0x1.5_5576_6647_2e04p-3)

        ## Step 3' (alternative)
        u = (f+f)/(y+F)
        v = u*u
        q = u*v*@horner(v,
                        0x1.5_5555_5555_0286p-4,
                        0x1.9_99a0_bc71_2416p-7)

        ## Step 4
        return l_lead + (u + (q + l_trail))
    elseif x == 0.0
        -Inf
    else
        NaN
    end
end
