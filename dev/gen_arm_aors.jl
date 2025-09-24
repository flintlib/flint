#
#   Copyright (C) 2024 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

# Generating routines for r <- a OP b, where OP is either + or -.
#
# This generation was constructed with processors with Apple silicon in mind.
# Processors decoding less than 6 operations per cycle, or few store and load
# units may have worse performance.

r = "rp"
a = "ap"
b = "bp"
rp(ix::Int) = "[$r,#$ix*8]"
ap(ix::Int) = "[$a,#$ix*8]"
bp(ix::Int) = "[$b,#$ix*8]"

sx = "sx" # Return value for carry or borrow
CC = "CC"

sp = ["s$ix" for ix in 0:14] # Scrap registers

# Writes assembly that should be preprocessed by M4.
function aors(n::Int)
    _str = "PROLOGUE(flint_mpn_aors($n))\n"
    function ldr(s0::String, s1::String)
        _str *= "\tldr\t$s0, $s1\n"
    end
    function ldp(s0::String, s1::String, s2::String)
        _str *= "\tldp\t$s0, $s1, $s2\n"
    end
    function str(s0::String, s1::String)
        _str *= "\tstr\t$s0, $s1\n"
    end
    function stp(s0::String, s1::String, s2::String)
        _str *= "\tstp\t$s0, $s1, $s2\n"
    end
    function OP(s0::String, s1::String, s2::String)
        _str *= "\tOP\t$s0, $s1, $s2\n"
    end
    function OPC(s0::String, s1::String, s2::String)
        _str *= "\tOPC\t$s0, $s1, $s2\n"
    end
    function cset(s0::String, s1::String)
        _str *= "\tcset\t$s0, $s1\n"
    end

    sv = deepcopy(sp)
    s(ix::Int) = sv[ix + 1]
    function shift(sv::Vector{String})
        sv[(end - 3):end], sv[1:(end - 4)] = sv[1:4], sv[5:end]
    end

    ldp(    s(0), s(2), ap(0))
    ldp(    s(1), s(3), bp(0))
    OP(     s(0), s(0), s(1))
    OPC(    s(2), s(2), s(3))
    stp(    s(0), s(2), rp(0))

    for ix in 1:(n ÷ 2 - 1)
        shift(sv)
        ldp(    s(0), s(2), ap(2 * ix))
        ldp(    s(1), s(3), bp(2 * ix))
        OPC(    s(0), s(0), s(1))
        OPC(    s(2), s(2), s(3))
        stp(    s(0), s(2), rp(2 * ix))
    end

    if n % 2 == 1
        ldr(    s(4), ap(n - 1))
        ldr(    s(5), bp(n - 1))
        OPC(    s(4), s(4), s(5))
        str(    s(4), rp(n - 1))
    end

    cset(   sx, CC)

    _str *= "\tret\nEPILOGUE()\n"

    return _str
end

function print_all_aors(nmax::Int = 16)
    for n in 2:nmax
        println(aors(n))
    end
end
