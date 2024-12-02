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

# Generating routines for r <- a OP 2^(cnt) * b, where OP is either + or -.

r = "rp"
a = "ap"
b = "bp"
cnt = "cnt"
rp(ix::Int) = "$ix*8($r)"
ap(ix::Int) = "$ix*8($a)"
bp(ix::Int) = "$ix*8($b)"

tnc = "tnc"
sx = "sx" # Return value for carry or borrow, i.e. %rax

R32(sx::String) = "R32($sx)"
R8(sx::String) = "R8($sx)"

s0 = "s0"
s1 = "s1"
s2 = "s2"
sp = ["s$ix" for ix in 0:2] # Scrap registers
s(ix::Int) = s[ix + 1]

# Writes assembly that should be preprocessed by M4.
function addrsh(n::Int)
    str = "\tALIGN(16)\nPROLOGUE(flint_mpn_addrsh_$n)\n"
    function mov(s0::String, s1::String)
        str *= "\tmov\t$s0, $s1\n"
    end
    function xor(s0::String, s1::String)
        str *= "\txor\t$s0, $s1\n"
    end
    function add(s0::String, s1::String)
        str *= "\tadd\t$s0, $s1\n"
    end
    function adc(s0::String, s1::String)
        str *= "\tadc\t$s0, $s1\n"
    end
    function sub(s0::String, s1::String)
        str *= "\tsub\t$s0, $s1\n"
    end
    function sbb(s0::String, s1::String)
        str *= "\tsbb\t$s0, $s1\n"
    end
    function shrx(s0::String, s1::String, s2::String)
        str *= "\tshrx\t$s0, $s1, $s2\n"
    end
    function shlx(s0::String, s1::String, s2::String)
        str *= "\tshlx\t$s0, $s1, $s2\n"
    end
    function lea(t::Tuple{String, String}, s1::String)
        str *= "\tlea\t($(t[1]), $(t[2])), $s1\n"
    end
    function setc(s0::String)
        str *= "\tsetc\t$s0\n"
    end

    # Initialize variables
    xor(    R32(tnc), R32(tnc))
    sub(    cnt, tnc)   # This is modulo 64, so -n = 64 - n.
    xor(    R32(sx), R32(sx))

    # f_a assumes s1 contains ix*8(bp)
    function f_a(ix::Int)
        if ix == 0
            shrx(   cnt, bp(0), s0)
            mov(    bp(ix + 1), s1)
        elseif ix == n - 1
            shrx(   cnt, s1, s0)
        else
            shrx(   cnt, s1, s0)
            mov(    bp(ix + 1), s1)
        end
    end # s0, s1 used
    function f_b(ix::Int)
        if ix != n - 1
            shlx(   tnc, s1, s2)
            lea(    (s0, s2), s2)
        end
    end # s1, s2 used
    function f_c(ix::Int)
        if ix == 0
            add(    ap(ix), s2)
            mov(    s2, rp(ix))
        elseif ix == n - 1
            adc(    ap(ix), s2)
            mov(    s2, rp(ix))
        else
            adc(    ap(ix), s0)
            mov(    s0, rp(ix))
        end
    end # nothing used

    # We interleave as follows:
    f_a(0)
    f_b(0)
    for ix in 1:(n - 1)
        f_a(ix + 0)
        f_c(ix - 1)
        f_b(ix + 0)
    end
    f_c(n - 1)

    setc(   R8(sx))

    str *= "\tret\nEPILOGUE()\n"

    return str
end

function print_all_addrsh(nmax::Int = 16)
    for n in 2:nmax
        println(addrsh(n))
    end
end
