
# A Julia translation of the log function of Yeppp:
#  https://bitbucket.org/MDukhan/yeppp/src/b8db687e912bbd7e2a26dd93c348ef2d7a5febdc/library/headers/yepBuiltin.h?at=default#cl-1128

# Available under following terms:
#
#                          Yeppp! library header
#
# This file is part of Yeppp! library and licensed under the New BSD license.
#
# Copyright (C) 2010-2012 Marat Dukhan
# Copyright (C) 2012-2013 Georgia Institute of Technology
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Georgia Institute of Technology nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function log_yeppp(x::Float64)
    xu = reinterpret(Uint64,x)
    if xu >= 0x0010_0000_0000_0000
        k = int(xu >> 52) & 0x07ff - 1023
        s = reinterpret(Float64,(xu & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000)
    else
        offset = leading_zeros(xu)-11
        k = -1022 - offset
        s = reinterpret(Float64,((xu << offset) & 0x000f_ffff_ffff_ffff) 
                        | 0x3ff0_0000_0000_0000)
    end
    if s >= 0x1.6A09E667F3BCDp+0 # sqrt(2)
        k += 1
        s *= 0.5
    end
    fk = Float64(k)

    t = s-1.0
    pt = @horner(t,
            -0x1.FFFFFFFFFFFF2p-2,
            0x1.5555555555103p-2,
            -0x1.00000000013C7p-2,
            0x1.9999999A43E4Fp-3,
            -0x1.55555554A6A2Bp-3,
            0x1.249248DAE4B2Ap-3,
            -0x1.FFFFFFBD8606Dp-4,
            0x1.C71C90DB06248p-4,
            -0x1.9999C5BE751E3p-4,
            0x1.745980F3FB889p-4,
            -0x1.554D5ACD502ABp-4,
            0x1.3B4ED39194B87p-4,
            -0x1.25480A82633AFp-4,
            0x1.0F23916A44515p-4,
            -0x1.EED2E2BB64B2Ep-5,
            0x1.EA17E14773369p-5,
            -0x1.1654764F478ECp-4,
            0x1.0266CD08DB2F2p-4,
            -0x1.CC4EC078138E3p-6)

    rf = t*(t*pt)+t
    f = fk*0x1.62E42FEFA3800p-1 + (fk*0x1.EF35793C76730p-45 + rf)
    if xu < 0
        NaN
    elseif xu == 0
        -Inf
    elseif !(x < Inf)
        x
    else
        f
    end
end

@vectorize_1arg Real log_yeppp
