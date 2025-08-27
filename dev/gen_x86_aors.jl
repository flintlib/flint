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
# This generation was constructed with processors with descent schedulers in
# mind.

r = "rp"
a = "ap"
b = "bp"
rp(ix::Int) = "$ix*8($r)"
ap(ix::Int) = "$ix*8($a)"
bp(ix::Int) = "$ix*8($b)"

sx = "sx" # Return value for carry or borrow, i.e. %rax

R32(sx::String) = "R32($sx)"
R8(sx::String) = "R8($sx)"

sp = ["s$ix" for ix in 0:4] # Scrap registers

# Writes assembly that should be preprocessed by M4.
function aors(n::Int)
    str = "\tALIGN(16)\nPROLOGUE(flint_mpn_aors($n))\n"
    function mov(s0::String, s1::String)
        str *= "\tmov\t$s0, $s1\n"
    end
    function xor(s0::String, s1::String)
        str *= "\txor\t$s0, $s1\n"
    end
    function OP(s0::String, s1::String)
        str *= "\tOP\t$s0, $s1\n"
    end
    function OPC(s0::String, s1::String)
        str *= "\tOPC\t$s0, $s1\n"
    end
    function setc(s0::String)
        str *= "\tsetc\t$s0\n"
    end

    sv = deepcopy(sp)
    s(ix::Int) = sv[ix + 1]
    function shift(sv::Vector{String})
        sv[end], sv[1:end - 1] = sv[1], sv[2:end]
    end

    mov(    ap(0), s(0))

    mov(    ap(1), s(1))
    xor(    R32(sx), R32(sx))
    OP(     bp(0), s(0))
    mov(    s(0), rp(0))

    for ix in 1:(n - 2)
        shift(sv)
        mov(    ap(ix + 1), s(1))
        OPC(    bp(ix), s(0))
        mov(    s(0), rp(ix))
    end

    OPC(    bp(n - 1), s(1))
    mov(    s(1), rp(n - 1))
    setc(   R8(sx))

    str *= "\tret\nEPILOGUE()\n"

    return str
end

function print_all_aors(nmax::Int = 16)
    for n in 2:nmax
        println(aors(n))
    end
end
