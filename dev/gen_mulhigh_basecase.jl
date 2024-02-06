#  AMD64 mpn_addmul_1 optimised for Intel Broadwell.
#
#  Copyright 2015, 2017 Free Software Foundation, Inc.
#
#  This file is part of the GNU MP Library.
#
#  The GNU MP Library is free software; you can redistribute it and/or modify
#  it under the terms of either:
#
#    * the GNU Lesser General Public License as published by the Free
#      Software Foundation; either version 3 of the License, or (at your
#      option) any later version.
#
#  or
#
#    * the GNU General Public License as published by the Free Software
#      Foundation; either version 2 of the License, or (at your option) any
#      later version.
#
#  or both in parallel, as here.
#
#  The GNU MP Library is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  You should have received copies of the GNU General Public License and the
#  GNU Lesser General Public License along with the GNU MP Library.  If not,
#  see https://www.gnu.org/licenses/.
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

# NOTE: GMP's addmul_1 code is used in mulhigh_basecase.

_regs = ["%rdi", "%rsi", "%rdx", "%rcx", "%r8", "%r9", "%r10", "%r11", "%rax"]
__regs = ["%rbx", "%rbp", "%r12", "%r13", "%r14", "%r15"]

function R8(reg::String)
    if reg == "%rax"
        return "%al"
    elseif reg == "%rbx"
        return "%bl"
    elseif reg == "%rcx"
        return "%cl"
    elseif reg == "%rdx"
        return "%dl"
    elseif reg == "%rsp"
        return "%spl"
    elseif reg == "%rbp"
        return "%bpl"
    elseif reg == "%rsi"
        return "%sil"
    elseif reg == "%rdi"
        return "%dil"
    elseif reg == "%r8"
        return "%r8b"
    elseif reg == "%r9"
        return "%r9b"
    elseif reg == "%r10"
        return "%r10b"
    elseif reg == "%r11"
        return "%r11b"
    elseif reg == "%r12"
        return "%r12b"
    elseif reg == "%r13"
        return "%r13b"
    elseif reg == "%r14"
        return "%r14b"
    elseif reg == "%r15"
        return "%r15b"
    else
        return "hejhoppgummi"
    end
end

function R32(reg::String)
    if reg == "%rax"
        return "%eax"
    elseif reg == "%rbx"
        return "%ebx"
    elseif reg == "%rcx"
        return "%ecx"
    elseif reg == "%rdx"
        return "%edx"
    elseif reg == "%rsp"
        return "%esp"
    elseif reg == "%rbp"
        return "%ebp"
    elseif reg == "%rsi"
        return "%esi"
    elseif reg == "%rdi"
        return "%edi"
    elseif reg == "%r8"
        return "%r8d"
    elseif reg == "%r9"
        return "%r9d"
    elseif reg == "%r10"
        return "%r10d"
    elseif reg == "%r11"
        return "%r11d"
    elseif reg == "%r12"
        return "%r12d"
    elseif reg == "%r13"
        return "%r13d"
    elseif reg == "%r14"
        return "%r14d"
    elseif reg == "%r15"
        return "%r15d"
    else
        return "hejhoppgummi"
    end
end

###############################################################################
# Preamble
###############################################################################

GMPaddmul1_copyright = "#  AMD64 mpn_addmul_1 optimised for Intel Broadwell.
#
#  Copyright 2015, 2017 Free Software Foundation, Inc.
#
#  This file is part of the GNU MP Library.
#
#  The GNU MP Library is free software; you can redistribute it and/or modify
#  it under the terms of either:
#
#    * the GNU Lesser General Public License as published by the Free
#      Software Foundation; either version 3 of the License, or (at your
#      option) any later version.
#
#  or
#
#    * the GNU General Public License as published by the Free Software
#      Foundation; either version 2 of the License, or (at your option) any
#      later version.
#
#  or both in parallel, as here.
#
#  The GNU MP Library is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  You should have received copies of the GNU General Public License and the
#  GNU Lesser General Public License along with the GNU MP Library.  If not,
#  see https://www.gnu.org/licenses/.
#
#   Copyright (C) 2024 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#\n"

copyright = "#
#   Copyright (C) 2024 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#\n"

preamble = "include(`config.m4')dnl\ndnl\n.text\n"

function function_pre_post(funname::String)
    pre = ".global\tFUNC($funname)
.p2align\t5, 0x90
TYPE($funname)

FUNC($funname):
\t.cfi_startproc\n"

    post = ".$(funname)_end:
SIZE($(funname), .$(funname)_end)
.cfi_endproc\n"

    return (pre, post)
end

###############################################################################
# mulhigh, hardcoded
###############################################################################

function mulhigh_1()
    r0 = _regs[9] # Important that r0 is rax
    r1 = _regs[4]

    res = _regs[1]
    ap = _regs[2]
    bp = _regs[3]

    body = ""

    body *= "\tmov\t0*8($bp), %rdx\n"
    body *= "\tmulx\t0*8($ap), $r0, $r1\n"
    body *= "\tmov\t$r1, 0*8($res)\n"

    return body * "\n\tret\n"
end

function mulhigh_2()
    r0 = _regs[9]
    r1 = _regs[4]
    r2 = _regs[5]

    sc = _regs[6]
    zr = _regs[7]

    res = _regs[1]
    ap = _regs[2]
    bp_old = _regs[3]
    b1 = _regs[8]

    body = ""

    body *= "\tmov\t1*8($bp_old), $b1\n"
    body *= "\tmov\t0*8($bp_old), %rdx\n"
    body *= "\txor\t$(R32(zr)), $(R32(zr))\n"
    body *= "\tmulx\t0*8($ap), $sc, $r0\n"
    body *= "\tmulx\t1*8($ap), $sc, $r1\n"
    body *= "\tadcx\t$sc, $r0\n"
    body *= "\tadcx\t$zr, $r1\n"

    body *= "\tmov\t$b1, %rdx\n"
    body *= "\tmulx\t0*8($ap), $sc, $r2\n"
    body *= "\tadcx\t$sc, $r0\n"
    body *= "\tadcx\t$r2, $r1\n"
    body *= "\tmulx\t1*8($ap), $sc, $r2\n"
    body *= "\tadox\t$sc, $r1\n"
    body *= "\tadox\t$zr, $r2\n"
    body *= "\tadcx\t$zr, $r2\n"
    body *= "\tmov\t$r1, 0*8($res)\n"
    body *= "\tmov\t$r2, 1*8($res)\n"

    return body * "\n\tret\n"
end

# When n = 9, push res to stack.
# When n = 10, push res and rsp to xmm register.
# When n = 11, do it like n = 10, and also use bp as r(11) as r(11) is only
# really needed in the last step.
# When n = 12, do it like n = 11 but do not use a zero register.
function mulhigh(n::Int; debug::Bool = false)
    if n < 1
        error()
    elseif n == 1
        return mulhigh_1()
    elseif n == 2
        return mulhigh_2()
    elseif n <= 12
        # Continue
    else
        error()
    end

    if debug
        res = "res"
        ap = "ap"
        bp_old = "bp_old"
        bp = "bp"
        sc = "sc"
        zr = "zr"
    else
        res = _regs[1]
        ap = _regs[2]
        bp_old = _regs[3] # rdx
        bp = _regs[4] # rdx is used by mulx, so we need to switch register for bp

        if n != 12
            sc = _regs[5] # scrap register
            zr = (n < 10) ? _regs[6] : "%rsp" # zero
        else
            sc = "%rsp"
            zr = sc
        end

        if n < 9
            _r = [_regs[9]; _regs[7:8]; __regs[1:n - 2]]
        elseif n == 9
            _r = [_regs[9]; _regs[7:8]; __regs[1:6]; res]
        elseif n == 10
            _r = [_regs[9]; _regs[6:8]; __regs[1:6]; res]
        elseif n == 11
            _r = [_regs[9]; _regs[6:8]; __regs[1:6]; res; bp]
        elseif n == 12
            _r = [_regs[9]; _regs[5:8]; __regs[1:6]; res; bp]
        end
    end

    r(ix::Int) = debug ? "r$ix" : _r[ix + 1]

    # Push
    body = ""
    for ix in 1:min(n - 2, 6)
        body *= "\tpush\t$(__regs[ix])\n"
    end
    if n == 9
        body *= "\tpush\t$res\n"
    elseif n >= 10
        body *= "\tvmovq\t%rsp, %xmm0\n"
        body *= "\tvmovq\t$res, %xmm1\n"
    end
    body *= "\n"

    # Prepare
    body *= "\tmov\t$bp_old, $bp\n"
    body *= "\tmov\t0*8($bp_old), %rdx\n"
    if n != 12
        body *= "\txor\t$(R32(zr)), $(R32(zr))\n"
    end
    body *= "\n"

    # First multiplication chain
    body *= "\tmulx\t$(n - 2)*8($ap), $sc, $(r(0))\n"
    body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(1))\n"
    if n != 12
        body *= "\tadcx\t$sc, $(r(0))\n"
        body *= "\tadcx\t$zr, $(r(1))\n"
    else
        body *= "\tadd\t$sc, $(r(0))\n"
        body *= "\tadc\t\$0, $(r(1))\n"
        body *= "\ttest\t%al, %al\n"
    end
    body *= "\n"

    # Intermediate multiplication chains
    for ix in 1:min(n - 2, (n != 12) ? 8 : 9)
        body *= "\tmov\t$ix*8($bp), %rdx\n"

        body *= "\tmulx\t$(n - 2 - ix)*8($ap), $sc, $(r(ix + 2))\n"
        body *= "\tmulx\t$(n - 1 - ix)*8($ap), $sc, $(r(ix + 1))\n"
        body *= "\tadcx\t$(r(ix + 2)), $(r(0))\n"
        body *= "\tadox\t$sc, $(r(0))\n"
        body *= "\tadcx\t$(r(ix + 1)), $(r(1))\n"

        for jx in 1:ix - 1
            body *= "\tmulx\t$(n - 1 - ix + jx)*8($ap), $sc, $(r(ix + 1))\n"
            body *= "\tadox\t$sc, $(r(jx + 0))\n"
            body *= "\tadcx\t$(r(ix + 1)), $(r(jx + 1))\n"
        end

        body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(ix + 1))\n"
        body *= "\tadox\t$sc, $(r(ix + 0))\n"
        if n == 12
            body *= "\tmov\t\$0, $(R32(zr))\n"
        end
        body *= "\tadcx\t$zr, $(r(ix + 1))\n"
        body *= "\tadox\t$zr, $(r(ix + 1))\n"

        body *= "\n"
    end

    if n >= 11
        N = n - 1
        body *= "\tmov\t$(n - 2)*8($bp), %rdx\n"
        for ix in 0:n - 2
            body *= "\tmulx\t$ix*8($ap), $sc, $(r(N))\n"
            if ix == 0
                body *= "\tadcx\t$(r(N)), $(r(ix + 0))\n"
            else
                body *= "\tadox\t$sc, $(r(ix - 1))\n"
                body *= "\tadcx\t$(r(N)), $(r(ix + 0))\n"
            end
        end
        body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(N))\n"
        body *= "\tadox\t$sc, $(r(N - 1))\n"
        if n == 12
            body *= "\tmov\t\$0, $(R32(zr))\n"
        end
        body *= "\tadcx\t$zr, $(r(N))\n"
        body *= "\tadox\t$zr, $(r(N))\n"
        body *= "\n"
    end

    # Last multiplication chain
    body *= "\tmov\t$(n - 1)*8($bp), %rdx\n"
    for ix in 0:n - 2
        body *= "\tmulx\t$ix*8($ap), $sc, $(r(n))\n"
        if ix % 2 == 0
            body *= "\tadcx\t$sc, $(r(ix + 0))\n"
            body *= "\tadcx\t$(r(n)), $(r(ix + 1))\n"
        else
            body *= "\tadox\t$sc, $(r(ix + 0))\n"
            body *= "\tadox\t$(r(n)), $(r(ix + 1))\n"
        end
    end
    body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(n))\n"
    if (n - 1) % 2 == 0
        body *= "\tadcx\t$sc, $(r(n - 1))\n"
    else
        body *= "\tadox\t$sc, $(r(n - 1))\n"
    end
    if n == 12
        body *= "\tmov\t\$0, $(R32(zr))\n"
    end
    body *= "\tadcx\t$zr, $(r(n))\n"
    if n == 9
        # Use scrap register for storing pointer to res
        res = sc
        body *= "\tpop\t$res\n"
    end
    body *= "\tadox\t$zr, $(r(n))\n"
    body *= "\n"

    if n == 10 || n == 11
        res, zr = sc, "error zr"
        body *= "\tvmovq\t%xmm1, $res\n"
        body *= "\tvmovq\t%xmm0, %rsp\n"
    elseif n == 12
        res = sc
        body *= "\tvmovq\t%xmm1, $res\n"
    end

    # Store result
    for ix in 1:n
        body *= "\tmov\t$(r(ix)), $(ix - 1)*8($res)\n"
    end
    body *= "\n"

    # Pop
    if n == 12
        body *= "\tvmovq\t%xmm0, %rsp\n"
    end
    for ix in min(n - 2, 6):-1:1
        body *= "\tpop\t$(__regs[ix])\n"
    end
    body *= "\n"

    if debug
        print(body * "\tret\n")
    else
        return body * "\tret\n"
    end
end

###############################################################################
# Initial triangle
###############################################################################
# ap + ap_os ~= ap + n - 1
# bp + bp_os ~= bp
# Triangle with side of size k.
# Assumes that bp[0] is already in rdx.
# r0 should be rax.
function macro_triangle(k::Int)
    pre = ".macro\ttr_$(k) ap=$(_regs[2]), ap_os=0, bp=$(_regs[5]), bp_os=0, "
    for ix in 0:k
        pre *= "r$ix, "
    end
    pre *= "sc, zr, zr32\n"
    post = ".endm\n"

    body = ""
    body *= "\tmov\t(\\bp_os + 0)*8(\\bp), %rdx\n"
    body *= "\tmulx\t(\\ap_os - 1)*8(\\ap), \\r0, \\r0\n"
    body *= "\tmulx\t(\\ap_os - 0)*8(\\ap), \\sc, \\r1\n"
    body *= "\tadd\t\\sc, \\r0\n"
    body *= "\tadc\t\$0, \\r1\n"
    body *= "\txor\t\\zr32, \\zr32\n"
    body *= "\n"

    for ix in 1:k - 1
        body *= "\tmov\t(\\bp_os + $ix)*8(\\bp), %rdx\n"
        body *= "\tmulx\t(\\ap_os - $(ix + 1))*8(\\ap), \\sc, \\sc\n"
        body *= "\tadcx\t\\sc, \\r0\n"
        for jx in 1:ix
            body *= "\tmulx\t(\\ap_os - $(ix - jx + 1))*8(\\ap), \\sc, \\r$(ix + 1)\n"
            body *= "\tadox\t\\sc, \\r$(jx - 1)\n"
            body *= "\tadcx\t\\r$(ix + 1), \\r$(jx - 0)\n"
        end
        body *= "\tmulx\t(\\ap_os - 0)*8(\\ap), \\sc, \\r$(ix + 1)\n"
        body *= "\tadox\t\\sc, \\r$(ix + 0)\n"
        body *= "\tadcx\t\\zr, \\r$(ix + 1)\n"
        body *= "\tadox\t\\zr, \\r$(ix + 1)\n"
        if ix != k - 1
            body *= "\n"
        end
    end

    return pre * body * post
end

# Scheme for n = 11:
#
# a/b 0  1  2  3  4  5  6  7  8  9 10
#   +--------------------------------
#  0|                            h  x
#  1|                         h  x  x
#  2|                      h  x  x  x
#  3|                   h  x  x  x  x
#  4|                h  x  x  x  x  x
#  5|          ^  h  x  x  x  x  x  x
#  6|          h  x  x  x  x  x  x  x
#  7|       h  x  ø  x  x  x  x  x  x <- (ap, bp) will initially point to ø
#  8|    h  x  x  x  x  x  x  x  x  x
#  9| h  x  x  x  x  x  x  x  x  x  x
# 10| x  x  x  x  x  x  x  x  x  x  x
#
# Scheme for n = 6:
#
# a/b 0  1  2  3  4  5
#   +-----------------
#  0|             h  x
#  1|          h  x  x
#  2|       h  x  ø  x <- (ap, bp) will initially point to ø
#  3|    h  x  x  x  x
#  4| h  x  x  x  x  x
#  5| x  x  x  x  x  x
#
# Scheme for n = 5:
#
# a/b 0  1  2  3  4
#   +--------------
#  0|          h  x
#  1|       h  x  ø <- (ap, bp) will initially point to ø
#  2|    h  x  x  x
#  3| h  x  x  x  x
#  4| x  x  x  x  x
#
# NOTE: The addmul chains in this file uses code from GMP's addmul_1 for
# Broadwell.
#
# NOTE: Assumes that n > k + 1.
function mulhigh_basecase(k::Int)
    # Needed registers with an initial triangle of size k:
    # - rp
    # - ap
    # - bp
    # - n
    # - rdx for mulx
    # - rax for ret
    #
    # - For initial triangle:
    #   * k registers
    #   * sc
    #   * zr
    #
    # - For memory addmul-chains:
    #   * 4 registers
    #   * m
    #   * m_save
    #
    # Hence, we need at least 6 + max(k + 2, 6) registers. We set k = 4, so we
    # need 12 registers.
    if k != 4
        error()
    end
    rp = _regs[1]
    ap = _regs[2]
    bp_old = _regs[3] # rdx
    n_old = _regs[4] # rcx
    bp = _regs[5] # rdx is used by mulx
    ret = _regs[9] # rax for storing the lowest limb

    body = "\tmov\t$bp_old, $bp\n"
    body *= "\tlea\t-$(k)*8($ap,$n_old,8), $ap\n" # ap += n - k
    body *= "\tlea\t$(k)*8($bp), $bp\n" # bp += k

    ###########################################################################
    # Push
    ###########################################################################
    for reg in __regs[1:3]
        body *= "\tpush\t$reg\n"
    end
    body *= "\n"

    ###########################################################################
    # Do initial triangle
    ###########################################################################
    r = [_regs[8]; __regs[1:3]]
    sc = _regs[6] # scrap register
    zr = _regs[7] # zero

    # Do triangle
    ap_os = k - 1
    bp_os = -k
    body *= "\ttr_$(k)\t$ap, $ap_os, $bp, $bp_os, $ret, "
    for ix in 1:k
        body *= "$(r[ix]), "
    end
    body *= "$sc, $zr, $(R32(zr))\n"

    # Push into rp
    for ix in 1:k
        body *= "\tmov\t$(r[ix]), $(ix - 1)*8($rp)\n"
    end
    body *= "\n"

    # Declare r, sc and zr dead
    r, sc, zr = "dead r", "dead sc", "dead zr"

    ###########################################################################
    # Addmul chains
    ###########################################################################
    # Set n
    n = _regs[6]
    body *= "\tlea\t-1($n_old), $n\n"

    # Set m and m_save
    m, m_save = _regs[4], _regs[7]
    body *= "\tmov\t\$$(k), $(R32(m_save))\n"
    body *= "\tmov\t\$$(k), $(R32(m))\n"

    # Declare four scrap registers
    s0, s1, s2, s3 = __regs[1], __regs[2], __regs[3], _regs[8]

    # Loop
    body *= "\t.align\t32, 0x90\n"
    body *= ".Lloop:"
    body *= "\tmov\t0*8($bp), %rdx\n"

    # Start computing stuff contributing to ret
    body *= "\tmulx\t-2*8($ap), $s1, $s1\n"

    # Set s0 to m % 8 and set m to m / 8
    body *= ".Lfin:"
    body *= "\tmov\t$(R32(m)), $(R32(s0))\n"
    body *= "\tshr\t\$3, $m\n"
    body *= "\tand\t\$7, $(R32(s0))\n" # Also resets flags from shr

    # Continue computing stuff contributing to ret
    body *= "\tmulx\t-1*8($ap), $s2, $s3\n"
    body *= "\tadcx\t$s1, $ret\n"
    body *= "\tadox\t$s2, $ret\n"
    body *= "\tmov\t$s3, $s1\n" # Either s1 or s3 is used before entering loop

    # Set s2 to pointer of .Ljmptab
    body *= "\tlea\t.Ljmptab(%rip), $s2\n"

    # Jump to label accordingly
    body *= "ifdef(`PIC',
`\tmovslq\t($s2,$s0,4), $s0
\tlea\t($s0,$s2), $s2
\tjmp\t*$s2
',`
\tjmp\t*($s2,$s0,8)
')\n"

    body *= "\t.section\t.data.rel.ro.local,\"a\",@progbits
\t.align\t8, 0x90
ifdef(`PIC',
`.Ljmptab:
\t.long\t.Lp0-.Ljmptab
\t.long\t.Lp1-.Ljmptab
\t.long\t.Lp2-.Ljmptab
\t.long\t.Lp3-.Ljmptab
\t.long\t.Lp4-.Ljmptab
\t.long\t.Lp5-.Ljmptab
\t.long\t.Lp6-.Ljmptab
\t.long\t.Lp7-.Ljmptab',
`.Ljmptab:
\t.quad\t.Lp0
\t.quad\t.Lp1
\t.quad\t.Lp2
\t.quad\t.Lp3
\t.quad\t.Lp4
\t.quad\t.Lp5
\t.quad\t.Lp6
\t.quad\t.Lp7')
\t.text\n"
    body *= "\n"

    # Precompute first elements to enter addmul loop
    body *= ".Lp0:"
    body *= "\tmulx\t0*8($ap), $s0, $s1\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= "\tlea\t-1*8($ap), $ap\n"
    body *= "\tlea\t-1*8($ap), $rp\n"
    body *= "\tlea\t-1($m), $m\n"
    body *= "\tjmp\t.Lam0\n"
    body *= ".Lp1:"
    body *= "\tmulx\t0*8($ap), $s2, $s3\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tjmp\t.Lam1\n"
    body *= ".Lp2:"
    body *= "\tmulx\t0*8($ap), $s0, $s1\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= "\tlea\t1*8($ap), $ap\n"
    body *= "\tlea\t1*8($rp), $rp\n"
    body *= "\tjmp\t.Lam2\n"
    body *= ".Lp3:"
    body *= "\tmulx\t0*8($ap), $s2, $s3\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tlea\t2*8($ap), $ap\n"
    body *= "\tlea\t-6*8($rp), $rp\n"
    body *= "\tjmp\t.Lam3\n"
    body *= ".Lp4:"
    body *= "\tmulx\t0*8($ap), $s0, $s1\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= "\tlea\t3*8($ap), $ap\n"
    body *= "\tlea\t-5*8($rp), $rp\n"
    body *= "\tjmp\t.Lam4\n"
    # For .Lp5, see below
    body *= ".Lp6:"
    body *= "\tmulx\t0*8($ap), $s0, $s1\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= "\tlea\t5*8($ap), $ap\n"
    body *= "\tlea\t-3*8($rp), $rp\n"
    body *= "\tjmp\t.Lam6\n"
    body *= ".Lp7:"
    body *= "\tmulx\t0*8($ap), $s2, $s3\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tlea\t-2*8($ap), $ap\n"
    body *= "\tlea\t-2*8($rp), $rp\n"
    body *= "\tjmp\t.Lam7\n"
    body *= "\n"
    body *= ".Lp5:"
    body *= "\tmulx\t0*8($ap), $s2, $s3\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tlea\t4*8($ap), $ap\n"
    body *= "\tlea\t-4*8($rp), $rp\n"
    # body *= "\tjmp\t.Lam5\n"

    body *= "\t.align\t32, 0x90\n"
    body *= ".Lam5:"
    body *= "\tmulx\t-3*8($ap), $s0, $s1\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= "\tadox\t4*8($rp), $s2\n"
    body *= "\tmov\t$s2, 4*8($rp)\n"
    body *= ".Lam4:"
    body *= "\tmulx\t-2*8($ap), $s2, $s3\n"
    body *= "\tadox\t5*8($rp), $s0\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tmov\t$s0, 5*8($rp)\n"
    body *= ".Lam3:"
    body *= "\tadox\t6*8($rp), $s2\n"
    body *= "\tmulx\t-1*8($ap), $s0, $s1\n"
    body *= "\tmov\t$s2, 6*8($rp)\n"
    body *= "\tlea\t8*8($rp), $rp\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= ".Lam2:"
    body *= "\tmulx\t0*8($ap), $s2, $s3\n"
    body *= "\tadox\t-1*8($rp), $s0\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tmov\t$s0, -1*8($rp)\n"
    if m != "%rcx"
        error()
    end
    body *= "\tjrcxz\t.Lend\n"
    body *= ".Lam1:"
    body *= "\tmulx\t1*8($ap), $s0, $s1\n"
    body *= "\tadox\t0*8($rp), $s2\n"
    body *= "\tlea\t-1($m), $m\n"
    body *= "\tmov\t$s2, 0*8($rp)\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= ".Lam0:"
    body *= "\tmulx\t2*8($ap), $s2, $s3\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tadox\t1*8($rp), $s0\n"
    body *= "\tmov\t$s0, 1*8($rp)\n"
    body *= ".Lam7:"
    body *= "\tmulx\t3*8($ap), $s0, $s1\n"
    body *= "\tlea\t8*8($ap), $ap\n"
    body *= "\tadcx\t$s3, $s0\n"
    body *= "\tadox\t2*8($rp), $s2\n"
    body *= "\tmov\t$s2, 2*8($rp)\n"
    body *= ".Lam6:"
    body *= "\tmulx\t-4*8($ap), $s2, $s3\n"
    body *= "\tadox\t3*8($rp), $s0\n"
    body *= "\tadcx\t$s1, $s2\n"
    body *= "\tmov\t$s0, 3*8($rp)\n"
    body *= "\tjmp\t.Lam5\n"
    body *= "\n"

    # Post addmul-loop
    body *= ".Lend:"
    body *= "\tadox\t0*8($rp), $s2\n"
    body *= "\tadcx\t$m, $s3\n" # Assumes m = 0
    body *= "\tadox\t$m, $s3\n"
    body *= "\tlea\t1($m_save), $m\n" # Increase m
    body *= "\tneg\t$m_save\n"
    body *= "\tlea\t1*8($bp), $bp\n" # Increase bp
    body *= "\tmov\t$s2, 0*8($rp)\n"
    body *= "\tmov\t$s3, 1*8($rp)\n"
    body *= "\tlea\t1*8($rp,$m_save,8), $rp\n" # Reset rp
    body *= "\tlea\t($ap,$m_save,8), $ap\n" # Reset ap
    body *= "\tneg\t$m_save\n"
    body *= "\tcmp\t$n, $m\n"
    body *= "\tlea\t1($m_save), $m_save\n" # Increase m_save

    # If m < n, continue looping
    body *= "\tjb\t.Lloop\n"

    # If m > n, exit
    body *= "\tja\t.Lexit\n"

    # Else, m == n and so we perform the last addmul chain
    body *= "\tmov\t0*8($bp), %rdx\n"
    body *= "\txor\t$(R32(s1)), $(R32(s1))\n" # Set s1 to zero
    body *= "\tjmp\t.Lfin\n"
    body *= "\n"

    ###########################################################################
    # Pop
    ###########################################################################
    body *= ".Lexit:"
    for reg in reverse(__regs[1:3])
        body *= "\tpop\t$reg\n"
    end
    body *= "\n"

    ###########################################################################
    # Return
    ###########################################################################
    body *= "\tret\n"
    return body
end

###############################################################################
# mulhigh, normalised
###############################################################################

function mulhigh_normalised_1()
    r0 = _regs[9] # Important that r0 is rax
    r1 = _regs[4]

    res = _regs[1]
    ap = _regs[2]
    bp = _regs[3]

    body = ""

    body *= "\tmov\t0*8($bp), %rdx\n"
    body *= "\tmulx\t0*8($ap), $r0, $r1\n"

    # Check if normalised
    body *= "\tmov\t\$0, %rdx\n"
    body *= "\ttest\t$r1, $r1\n"
    body *= "\tsetns\t$(R8("%rdx"))\n"
    body *= "\tjs\t.Lcontinue\n"

    # If not normalised, shift by one
    body *= "\tadd\t$r0, $r0\n"
    body *= "\tadc\t$r1, $r1\n"

    body *= ".Lcontinue:\n"
    body *= "\tmov\t$r1, 0*8($res)\n"

    return body * "\n\tret\n"
end

function mulhigh_normalised_2()
    r0 = _regs[9]
    r1 = _regs[4]
    r2 = _regs[5]

    sc = _regs[6]
    zr = _regs[7]

    res = _regs[1]
    ap = _regs[2]
    bp_old = _regs[3]
    b1 = _regs[8]

    body = ""

    body *= "\tmov\t1*8($bp_old), $b1\n"
    body *= "\tmov\t0*8($bp_old), %rdx\n"
    body *= "\txor\t$(R32(zr)), $(R32(zr))\n"
    body *= "\tmulx\t0*8($ap), $sc, $r0\n"
    body *= "\tmulx\t1*8($ap), $sc, $r1\n"
    body *= "\tadcx\t$sc, $r0\n"
    body *= "\tadcx\t$zr, $r1\n"

    body *= "\tmov\t$b1, %rdx\n"
    body *= "\tmulx\t0*8($ap), $sc, $r2\n"
    body *= "\tadcx\t$sc, $r0\n"
    body *= "\tadcx\t$r2, $r1\n"
    body *= "\tmulx\t1*8($ap), $sc, $r2\n"
    body *= "\tadox\t$sc, $r1\n"
    body *= "\tadox\t$zr, $r2\n"
    body *= "\tadcx\t$zr, $r2\n"

    # Check if normalised
    body *= "\tmov\t\$0, %rdx\n"
    body *= "\ttest\t$r2, $r2\n"
    body *= "\tsetns\t$(R8("%rdx"))\n"
    body *= "\tjs\t.Lcontinue\n"

    # If not normalised, shift by one
    body *= "\tadd\t$r0, $r0\n"
    body *= "\tadc\t$r1, $r1\n"
    body *= "\tadc\t$r2, $r2\n"

    body *= ".Lcontinue:\n"
    body *= "\tmov\t$r1, 0*8($res)\n"
    body *= "\tmov\t$r2, 1*8($res)\n"

    return body * "\n\tret\n"
end

function mulhigh_normalised(n::Int; debug::Bool = false)
    if n < 1
        error()
    elseif n == 1
        return mulhigh_normalised_1()
    elseif n == 2
        return mulhigh_normalised_2()
    elseif n <= 12
        # Continue
    else
        error()
    end

    if debug
        res = "res"
        ap = "ap"
        bp_old = "bp_old"
        bp = "bp"
        sc = "sc"
        zr = "zr"
    else
        res = _regs[1]
        ap = _regs[2]
        bp_old = _regs[3] # rdx
        bp = _regs[4] # rdx is used by mulx, so we need to switch register for bp

        if n != 12
            sc = _regs[5] # scrap register
            zr = (n < 10) ? _regs[6] : "%rsp" # zero
        else
            sc = "%rsp"
            zr = sc
        end

        if n < 9
            _r = [_regs[9]; _regs[7:8]; __regs[1:n - 2]]
        elseif n == 9
            _r = [_regs[9]; _regs[7:8]; __regs[1:6]; res]
        elseif n == 10
            _r = [_regs[9]; _regs[6:8]; __regs[1:6]; res]
        elseif n == 11
            _r = [_regs[9]; _regs[6:8]; __regs[1:6]; res; bp]
        elseif n == 12
            _r = [_regs[9]; _regs[5:8]; __regs[1:6]; res; bp]
        end
    end

    r(ix::Int) = debug ? "r$ix" : _r[ix + 1]

    # Push
    body = ""
    for ix in 1:min(n - 2, 6)
        body *= "\tpush\t$(__regs[ix])\n"
    end
    if n == 9
        body *= "\tpush\t$res\n"
    elseif n >= 10
        body *= "\tvmovq\t%rsp, %xmm0\n"
        body *= "\tvmovq\t$res, %xmm1\n"
    end
    body *= "\n"

    # Prepare
    body *= "\tmov\t$bp_old, $bp\n"
    body *= "\tmov\t0*8($bp_old), %rdx\n"
    if n != 12
        body *= "\txor\t$(R32(zr)), $(R32(zr))\n"
    end
    body *= "\n"

    # First multiplication chain
    body *= "\tmulx\t$(n - 2)*8($ap), $sc, $(r(0))\n"
    body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(1))\n"
    if n != 12
        body *= "\tadcx\t$sc, $(r(0))\n"
        body *= "\tadcx\t$zr, $(r(1))\n"
    else
        body *= "\tadd\t$sc, $(r(0))\n"
        body *= "\tadc\t\$0, $(r(1))\n"
        body *= "\ttest\t%al, %al\n"
    end
    body *= "\n"

    # Intermediate multiplication chains
    for ix in 1:min(n - 2, (n != 12) ? 8 : 9)
        body *= "\tmov\t$ix*8($bp), %rdx\n"

        body *= "\tmulx\t$(n - 2 - ix)*8($ap), $sc, $(r(ix + 2))\n"
        body *= "\tmulx\t$(n - 1 - ix)*8($ap), $sc, $(r(ix + 1))\n"
        body *= "\tadcx\t$(r(ix + 2)), $(r(0))\n"
        body *= "\tadox\t$sc, $(r(0))\n"
        body *= "\tadcx\t$(r(ix + 1)), $(r(1))\n"

        for jx in 1:ix - 1
            body *= "\tmulx\t$(n - 1 - ix + jx)*8($ap), $sc, $(r(ix + 1))\n"
            body *= "\tadox\t$sc, $(r(jx + 0))\n"
            body *= "\tadcx\t$(r(ix + 1)), $(r(jx + 1))\n"
        end

        body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(ix + 1))\n"
        body *= "\tadox\t$sc, $(r(ix + 0))\n"
        if n == 12
            body *= "\tmov\t\$0, $(R32(zr))\n"
        end
        body *= "\tadcx\t$zr, $(r(ix + 1))\n"
        body *= "\tadox\t$zr, $(r(ix + 1))\n"

        body *= "\n"
    end

    if n >= 11
        N = n - 1
        body *= "\tmov\t$(n - 2)*8($bp), %rdx\n"
        for ix in 0:n - 2
            body *= "\tmulx\t$ix*8($ap), $sc, $(r(N))\n"
            if ix == 0
                body *= "\tadcx\t$(r(N)), $(r(ix + 0))\n"
            else
                body *= "\tadox\t$sc, $(r(ix - 1))\n"
                body *= "\tadcx\t$(r(N)), $(r(ix + 0))\n"
            end
        end
        body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(N))\n"
        body *= "\tadox\t$sc, $(r(N - 1))\n"
        if n == 12
            body *= "\tmov\t\$0, $(R32(zr))\n"
        end
        body *= "\tadcx\t$zr, $(r(N))\n"
        body *= "\tadox\t$zr, $(r(N))\n"
        body *= "\n"
    end

    # Last multiplication chain
    body *= "\tmov\t$(n - 1)*8($bp), %rdx\n"
    for ix in 0:n - 2
        body *= "\tmulx\t$ix*8($ap), $sc, $(r(n))\n"
        if ix % 2 == 0
            body *= "\tadcx\t$sc, $(r(ix + 0))\n"
            body *= "\tadcx\t$(r(n)), $(r(ix + 1))\n"
        else
            body *= "\tadox\t$sc, $(r(ix + 0))\n"
            body *= "\tadox\t$(r(n)), $(r(ix + 1))\n"
        end
    end
    body *= "\tmulx\t$(n - 1)*8($ap), $sc, $(r(n))\n"
    if (n - 1) % 2 == 0
        body *= "\tadcx\t$sc, $(r(n - 1))\n"
    else
        body *= "\tadox\t$sc, $(r(n - 1))\n"
    end
    if n == 12
        body *= "\tmov\t\$0, $(R32(zr))\n"
    end
    body *= "\tadcx\t$zr, $(r(n))\n"
    if n == 9
        # Use scrap register for storing pointer to res
        res = sc
        body *= "\tpop\t$res\n"
    end
    body *= "\tadox\t$zr, $(r(n))\n"
    body *= "\n"

    # Check if normalised
    body *= "\tmov\t\$0, %rdx\n"
    body *= "\ttest\t$(r(n)), $(r(n))\n"
    body *= "\tsetns\t$(R8("%rdx"))\n"
    body *= "\tjs\t.Lcontinue\n"

    # If not normalised, shift by one
    body *= "\tadd\t$(r(0)), $(r(0))\n"
    for ix in 1:n
        body *= "\tadc\t$(r(ix)), $(r(ix))\n"
    end

    body *= ".Lcontinue:\n"
    if n == 10 || n == 11
        res, zr = sc, "error zr"
        body *= "\tvmovq\t%xmm1, $res\n"
        body *= "\tvmovq\t%xmm0, %rsp\n"
    elseif n == 12
        res = sc
        body *= "\tvmovq\t%xmm1, $res\n"
    end

    # Store result
    for ix in 1:n
        body *= "\tmov\t$(r(ix)), $(ix - 1)*8($res)\n"
    end
    body *= "\n"

    # Pop
    if n == 12
        body *= "\tvmovq\t%xmm0, %rsp\n"
    end
    for ix in min(n - 2, 6):-1:1
        body *= "\tpop\t$(__regs[ix])\n"
    end
    body *= "\n"

    if debug
        print(body * "\tret\n")
    else
        return body * "\tret\n"
    end
end

###############################################################################
# Generate file
###############################################################################

function gen_mulhigh_basecase(k::Int = 4, nofile::Bool = false)
    (pre, post) = function_pre_post("flint_mpn_mulhigh_n_basecase")
    macros = macro_triangle(k) * "\n"
    functionbody = mulhigh_basecase(k)

    str = "$GMPaddmul1_copyright\n$preamble\n$macros$pre$functionbody$post"

    if nofile
        print(str)
    else
        path = String(@__DIR__) * "/../src/mpn_extras/broadwell/mulhigh_n_basecase.asm"
        file = open(path, "w")
        write(file, str)
        close(file)
    end
end

function gen_mulhigh(m::Int, nofile::Bool = false)
    (pre, post) = function_pre_post("flint_mpn_mulhigh_$m")
    functionbody = mulhigh(m)

    str = "$copyright\n$preamble\n$pre$functionbody$post"

    if nofile
        print(str)
    else
        path = String(@__DIR__) * "/../src/mpn_extras/broadwell/mulhigh_$m.asm"
        file = open(path, "w")
        write(file, str)
        close(file)
    end
end

function gen_mulhigh_normalised(m::Int, nofile::Bool = false)
    (pre, post) = function_pre_post("flint_mpn_mulhigh_normalised_$m")
    functionbody = mulhigh_normalised(m)

    str = "$copyright\n$preamble\n$pre$functionbody$post"

    if nofile
        print(str)
    else
        path = String(@__DIR__) * "/../src/mpn_extras/broadwell/mulhigh_normalised_$m.asm"
        file = open(path, "w")
        write(file, str)
        close(file)
    end
end

function gen_all()
    gen_mulhigh_basecase()
    for m in 1:12
        gen_mulhigh(m)
    end
    for m in 1:12
        gen_mulhigh_normalised(m)
    end
end
