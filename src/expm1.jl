# computes 2^k exp(x) - 1 where |x| < ln(2)/2
# this is simply a translation of the last bit of openlib/src/s_expm1.c
# the only change is the handling of k == 0 to handle the extra precision term

function expm1_kernel(k::Int, x::Double64)
    hfx = 0.5*x.hi
    hxs = x.hi*hfx
    r1 = @horner(hxs,1.0,
                 -3.33333333333331316428e-02,
                 1.58730158725481460165e-03,
                 -7.93650757867487942473e-05,
                 4.00821782732936239552e-06,
                 -2.01099218183624371326e-07)
    
    t  = 3.0-r1*hfx    
    e  = hxs*((r1-t)/(6.0 - x.hi*t))
    if k == 0
        return x.hi - ((x.hi*(e-x.lo)-x.lo)-hxs) # modified to handle lo term.
    else
        twopk = ldexp(1.0,k)
	e  = (x.hi*(e-x.lo)-x.lo) - hxs
	if k == -1
            return 0.5*(x.hi-e)-0.5
        elseif k == 1
	    if x.hi < -0.25
                return -2.0*(e-(x.hi+0.5))
	    else 	      
                return 1.0+2.0*(x.hi-e)
	    end
	elseif k <= -2 || k > 56 # suffice to return exp(x)-1 
	    y = 1.0-(e-x.hi)
	    if k == 1024
                y = y*2.0*0x1p1023
	    else 
                y = y*twopk
	    end
	    return y-1.0
        end
	t = 1.0
	if k < 20 
            t = 1.0-ldexp(1.0,-k)
	    y = t-(e-x.hi)
	    y = y*twopk
	else
            t = ldexp(1.0,-k)
	    y = x.hi-(e+t)
	    y += 1.0
	    y = y*twopk
        end
    end
    y
end
