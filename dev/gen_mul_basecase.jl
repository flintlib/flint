#
#   Copyright (C) 2023 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 2.1 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

_regs = ["%rdi", "%rsi", "%rdx", "%rcx", "%r8", "%r9", "%r10", "%r11", "%rax"]

# We have omitted rsp as we cannot use it unless we push it to xmm register
__regs = ["%rbx", "%rbp", "%r12", "%r13", "%r14", "%r15"]

function reg_8_bit(reg::String)
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

function reg_32_bit(reg::String)
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

copyright = "#
#   Copyright (C) 2023 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 2.1 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#\n"

preamble = "include(`config.m4')dnl\ndnl\n.text\n"

function function_pre_post(funname::String)
    pre = ".global\tFUNC($funname)
.p2align\t4, 0x90
TYPE($funname)

FUNC($funname):
\t.cfi_startproc\n"

    post = ".$(funname)_end:
SIZE($(funname), .$(funname)_end)
.cfi_endproc\n"

    return (pre, post)
end

###############################################################################
# m = 1, n = 1
###############################################################################

function function_body_1(m::Int, n::Int = m)
    if m != 1 || n != 1
        error()
    end

    regs = ["%rcx", "%rax", "%r8", "%r9", "%r10", "%r11"]

    numregs = 2 # Number of registers used

    r0 = regs[1]
    r1 = regs[2] # Important that r1 is rax

    res_reg = "%rdi"
    ap_reg = "%rsi"
    bp_reg = "%rdx"

    body = ""

    body *= "\tmov\t($bp_reg), %rdx\n"
    body *= "\tmulx\t0*8($ap_reg), $r0, $r1\n"
    body *= "\tmov\t$r0, 0*8($res_reg)\n"
    body *= "\tmov\t$r1, 1*8($res_reg)\n"

    return body * "\n\tret\n"
end

###############################################################################
# m = 2, n = 1
###############################################################################

function function_body_2_1(m::Int, n::Int = 1)
    if m != 2 || n != 1
        error()
    end

    regs = ["%rcx", "%r8", "%rax", "%r9", "%r10", "%r11"]

    numregs = 5 # Number of registers used

    r0 = regs[1]
    r1 = regs[2]
    r2 = regs[3] # Important that r2 is rax
    r3 = regs[4]
    zero_reg = regs[numregs]
    zero_reg_32 = reg_32_bit(zero_reg)

    res_reg = "%rdi"
    ap_reg = "%rsi"
    bp_reg = "%rdx"

    body = ""

    body *= "\tmov\t0*8($bp_reg), %rdx\n"
    body *= "\txor\t$zero_reg_32, $zero_reg_32\n"
    body *= "\tmulx\t0*8($ap_reg), $r0, $r1\n"
    body *= "\tmulx\t1*8($ap_reg), $r3, $r2\n"
    body *= "\tadcx\t$r3, $r1\n"
    body *= "\tadcx\t$zero_reg, $r2\n"
    body *= "\tmov\t$r0, 0*8($res_reg)\n"
    body *= "\tmov\t$r1, 1*8($res_reg)\n"
    body *= "\tmov\t$r2, 2*8($res_reg)\n"

    return body * "\n\tret\n"
end

###############################################################################
# m = 2, n = 2
###############################################################################

function function_body_2(m::Int, n::Int = m)
    if m != 2 || n != 2
        error()
    end

    regs = ["%rcx", "%r8", "%r9", "%rax", "%r10", "%r11"]

    numregs = 6 # Number of registers used

    r0 = regs[1]
    r1 = regs[2]
    r2 = regs[3]
    r3 = regs[4] # Important that r3 is rax
    bp1_reg = regs[numregs - 1]
    zero_reg = regs[numregs]
    zero_reg_32 = reg_32_bit(zero_reg)

    res_reg = "%rdi"
    ap_reg = "%rsi"
    bp_reg = "%rdx"

    body = ""

    body *= "\tmov\t1*8($bp_reg), $bp1_reg\n"
    body *= "\tmov\t0*8($bp_reg), %rdx\n"
    body *= "\txor\t$zero_reg_32, $zero_reg_32\n"
    body *= "\tmulx\t0*8($ap_reg), $r0, $r1\n"
    body *= "\tmulx\t1*8($ap_reg), $r3, $r2\n"
    body *= "\tadcx\t$r3, $r1\n"
    body *= "\tmov\t$r0, 0*8($res_reg)\n"
    body *= "\tmov\t$bp1_reg, %rdx\n"
    body *= "\tmulx\t0*8($ap_reg), $r0, $r3\n"
    body *= "\tadox\t$r0, $r1\n"
    body *= "\tadcx\t$r3, $r2\n"
    body *= "\tmov\t$r1, 1*8($res_reg)\n"
    body *= "\tmulx\t1*8($ap_reg), $r0, $r3\n"
    body *= "\tadox\t$r0, $r2\n"
    body *= "\tadcx\t$zero_reg, $r3\n"
    body *= "\tadox\t$zero_reg, $r3\n"
    body *= "\tmov\t$r2, 2*8($res_reg)\n"
    body *= "\tmov\t$r3, 3*8($res_reg)\n"

    return body * "\n\tret\n"
end

###############################################################################
# mul_1 and addmul_1 macros
###############################################################################

# Assumes bp[0] is already in %rdx.
# This pushes directly into res, apart from the last limb which is left in r5.
function mul_1_macro(m::Int)
    if m < 3
        error()
    end

    pre = ".macro\tm$m res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, "
    pre *= "r0, r1, r2, r3, r4, r5\n"
    post = ".endm\n"

    r0 = "\\r0"
    r1 = "\\r1"
    r2 = "\\r2"
    r3 = "\\r3"
    r4 = "\\r4"
    r5 = "\\r5"

    # Make sure that last limb is pushed into r5
    if m % 2 == 0
        r4, r5 = r5, r4
    elseif m % 4 == 1
        r1, r5 = r5, r1
    end

    body = ""

    body *= "\tmulx\t(0+\\ap_offset)*8(\\ap), $r0, $r3\n"
    body *= "\tmulx\t(1+\\ap_offset)*8(\\ap), $r1, $r4\n"
    body *= "\tmulx\t(2+\\ap_offset)*8(\\ap), $r2, $r5\n"
    body *= "\tadd\t$r3, $r1\n"
    body *= "\tadc\t$r4, $r2\n"
    body *= "\tmov\t$r0, (0+\\res_offset)*8(\\res)\n"
    body *= "\tmov\t$r1, (1+\\res_offset)*8(\\res)\n"
    body *= "\tmov\t$r2, (2+\\res_offset)*8(\\res)\n"

    for ix in 1:(m - 1) ÷ 2 - 1
        i1 = 2 * ix + 1
        i2 = 2 * ix + 2
        body *= "\tmulx\t($i1+\\ap_offset)*8(\\ap), $r3, $r4\n"
        body *= "\tmulx\t($i2+\\ap_offset)*8(\\ap), $r0, $r1\n"
        body *= "\tadc\t$r3, $r5\n"
        body *= "\tadc\t$r4, $r0\n"
        body *= "\tmov\t$r5, ($i1+\\res_offset)*8(\\res)\n"
        body *= "\tmov\t$r0, ($i2+\\res_offset)*8(\\res)\n"
        (r0, r2), (r3, r4), (r1, r5) = (r2, r0), (r3, r4), (r5, r1)
    end

    if m % 2 == 1
        body *= "\tadc\t\$0, $r5\n"
    else
        body *= "\tmulx\t($(m - 1)+\\ap_offset)*8(\\ap), $r3, $r4\n"
        body *= "\tadc\t$r3, $r5\n"
        body *= "\tadc\t\$0, $r4\n"
        body *= "\tmov\t$r5, ($(m - 1)+\\res_offset)*8(\\res)\n"
    end

    return pre * body * post
end

# Assumes b is already in %rdx.
# Pushes least significant limb into res[res_offset].
# ip1 can be aliased with rN for N > 0 and scr2
# ip2 can be aliased with rN for N > 1
function mulM_macro(n::Int; chain::Bool = false)
    if n < 2 || n > 8
        error()
    end

    pre = ".macro\tm$(n)$(chain ? "_chain" : "") res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, "
    if chain
        pre *= "ip1, ip2, "
    end
    for jx in 0:n-1
        pre *= "r$jx, "
    end
    pre *= "scr1, scr2, zero\n"

    post = ".endm\n"

    body = ""
    body *= "\tmulx\t(0+\\ap_offset)*8(\\ap), \\scr1, \\r0\n"
    if chain
        body *= "\tadcx\t\\ip1, \\scr1\n"
        body *= "\tmov\t\\scr1, \\res_offset*8(\\res)\n"
    end
    body *= "\tmulx\t(1+\\ap_offset)*8(\\ap), \\scr2, \\r1\n"
    if !chain
        body *= "\tmov\t\\scr1, \\res_offset*8(\\res)\n"
    end
    body *= "\tadcx\t\\scr2, \\r0\n"
    if chain
        body *= "\tadox\t\\ip2, \\r0\n"
    end

    for jx in 2:n ÷ 2
        jxm = 2 * (jx - 1)
        body *= "\tmulx\t($(jxm)+\\ap_offset)*8(\\ap), \\scr1, \\r$jxm\n"
        body *= "\tmulx\t($(jxm + 1)+\\ap_offset)*8(\\ap), \\scr2, \\r$(jxm + 1)\n"
        body *= "\tadcx\t\\scr1, \\r$(jxm - 1)\n"
        body *= "\tadcx\t\\scr2, \\r$jxm\n"
    end

    if n % 2 == 1
        body *= "\tmulx\t($(n - 1)+\\ap_offset)*8(\\ap), \\scr1, \\r$(n - 1)\n"
        body *= "\tadcx\t\\scr1, \\r$(n - 2)\n"
    end

    body *= "\tadcx\t\\zero, \\r$(n - 1)\n"

    return pre * body * post
end

# Assumes b is already in %rdx.
# Pushes least significant limb into res[res_offset].
function addmulM_macro(n::Int)
    if n < 2 || n > 8
        error()
    end

    pre = ".macro\tam$(n) res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, "
    for jx in 0:n
        pre *= "r$jx, "
    end
    pre *= "scr, zero\n"
    post = ".endm\n"

    scr1 = "\\r$n"
    scr2 = "\\scr"
    if n % 2 == 1
        # If n is odd, we need to reorder scr1 and scr2 to make \scr occur at
        # the right place so that we do not need two scrap registers in the last
        # multiplication.
        tmp = scr1
        scr1 = scr2
        scr2 = tmp
    end

    body = ""
    body *= "\tmulx\t(0+\\ap_offset)*8(\\ap), $scr1, $scr2\n"
    body *= "\tadcx\t$scr1, \\r0\n"
    body *= "\tmov\t\\r0, \\res_offset*8(\\res)\n"

    for jx in 1:n - 1
        body *= "\tmulx\t($(jx)+\\ap_offset)*8(\\ap), \\r0, $(jx % 2 == 1 ? scr1 : scr2)\n"
        body *= "\tadcx\t$(jx % 2 == 0 ? scr1 : scr2), \\r$jx\n"
        body *= "\tadox\t\\r0, \\r$jx\n"
    end

    body *= "\tadcx\t\\zero, \\r$n\n"
    body *= "\tadox\t\\zero, \\r$n\n"

    return pre * body * post
end

###############################################################################
# m ≥ 3, n = 1
###############################################################################

function function_body_M_1(m::Int)
    if m < 3
        error()
    end

    body = "\tmov\t0*8(%rdx), %rdx\n"
    # Important that the last entry is rax to set return value
    body *= "\tm$m\t%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax\n"
    body *= "\tmov\t%rax, $m*8(%rdi)\n"

    return body * "\n\tret\n"
end

###############################################################################
# m > 8, n = 2
###############################################################################

function M_2_getN(m::Int)
    if m <= 8
        error()
    elseif m <= 9
        return 3
    elseif m <= 12
        return 4
    elseif m <= 15
        return 5
    elseif m <= 16
        return 6
    else
        error()
    end
end

# Needed macros: mN, amN_chain, mN_chain (if m ÷ N == 2),
# m(m % N)_chain (if m % N > 1) and am(m % N) (if m % N > 1).
function function_body_M_2(m::Int)
    if m <= 8
        error()
    end

    N = M_2_getN(m)

    # number of available registers is 10
    regs = ["%rax", "%r9", "%r10", "%r11", "%rbx", "%rbp", "%r12", "%r13", "%r14", "%r15"]
    preserve_index = 5 # Starting with %rbx, we have to push it to the stack

    used_registers = N + 3

    res = "%rdi"
    ap = "%rsi"
    bp = "%rdx"
    b0 = "%rcx"
    b1 = "%r8"

    if used_registers > length(regs)
        error()
    end

    pre = "\tmov\t0*8($bp), $b0\n"
    pre *= "\tmov\t1*8($bp), $b1\n"
    post = ""
    for jx in preserve_index:used_registers
        pre = pre * "\tpush\t$(regs[jx])\n"
        post = "\tpop\t$(regs[jx])\n" * post
    end

    rg = deepcopy(regs[1:used_registers])
    if m % N == 0 && (m ÷ N) % 2 == 1
        rg[1], rg[N + 1] = rg[N + 1], rg[1]
    elseif m % N == 1
        rg[1], rg[4] = rg[4], rg[1]
    elseif m % N > 1 && (m ÷ N) % 2 == 0
        rg[1], rg[N + 1] = rg[N + 1], rg[1]
    end

    scr = rg[N + 2]
    zero = rg[N + 3]

    body = ""

    body *= "\n"
    body *= "\tmov\t$b0, %rdx\n"
    body *= "\txor\t$(reg_32_bit(zero)), $(reg_32_bit(zero))\n"
    body *= "\tm$N\t$res, 0, $ap, 0, "
    for jx in 2:N
        body *= "$(rg[jx]), "
    end
    body *= "$(rg[1]), $(rg[N + 1]), $scr, $zero\n"

    body *= "\tmov\t$b1, %rdx\n"
    body *= "\tam$(N)\t$res, 1, $ap, 0, "
    for jx in 2:N
        body *= "$(rg[jx]), "
    end
    body *= "$(rg[1]), $(rg[N + 1]), $scr, $zero\n"
    for jx in 3:N
        body *= "\tmov\t$(rg[jx]), $(jx - 1)*8($res)\n"
    end

    for ix in 1:m ÷ N - 1
        body *= "\n"
        body *= "\tmov\t$b0, %rdx\n"
        body *= "\tm$(N)_chain\t$res, $(ix * N), $ap, $(ix * N), $(rg[1]), $(rg[N + 1]), "
        for jx in 2:N
            body *= "$(rg[jx]), "
        end
        body *= "$(rg[N + 1]), $scr, $(rg[1]), $zero\n"

        body *= "\tmov\t$b1, %rdx\n"
        body *= "\tam$(N)\t$res, $(ix * N + 1), $ap, $(ix * N), "
        for jx in 2:N
            body *= "$(rg[jx]), "
        end
        body *= "$(rg[N + 1]), $(rg[1]), $scr, $zero\n"
        for jx in 3:N
            body *= "\tmov\t$(rg[jx]), $(jx + ix * N - 1)*8($res)\n"
        end
        rg[1], rg[N + 1] = rg[N + 1], rg[1]
    end

    body *= "\n"
    if m % N == 0
        body *= "\tmov\t$(rg[1]), $(m)*8($res)\n"
        body *= "\tmov\t$(rg[N + 1]), $(m + 1)*8($res)\n"
    elseif m % N == 1
        body *= "\tmov\t$b0, %rdx\n"
        body *= "\tmulx\t$((m ÷ N) * N)*8($ap), $(rg[2]), $scr\n"
        body *= "\tmov\t$b1, %rdx\n"
        body *= "\tmulx\t$((m ÷ N) * N)*8($ap), $(rg[3]), $(rg[4])\n"
        body *= "\tadcx\t$(rg[1]), $(rg[2])\n"
        body *= "\tadcx\t$(rg[N + 1]), $(rg[3])\n"
        body *= "\tadcx\t$zero, $(rg[4])\n"
        body *= "\tadox\t$scr, $(rg[3])\n"
        body *= "\tadox\t$zero, $(rg[4])\n"
        body *= "\tmov\t$(rg[2]), $((m ÷ N) * N + 0)*8($res)\n"
        body *= "\tmov\t$(rg[3]), $((m ÷ N) * N + 1)*8($res)\n"
        body *= "\tmov\t$(rg[4]), $((m ÷ N) * N + 2)*8($res)\n"
    else
        body *= "\tmov\t$b0, %rdx\n"
        body *= "\tm$(m % N)_chain\t$res, $((m ÷ N) * N + 0), $ap, $((m ÷ N) * N), $(rg[1]), $(rg[N + 1]), "
        for jx in 2:m % N + 1
            body *= "$(rg[jx]), "
        end
        body *= "$scr, $(rg[1]), $zero\n"

        body *= "\tmov\t$b1, %rdx\n"
        body *= "\tam$(m % N)\t$res, $((m ÷ N) * N + 1), $ap, $((m ÷ N) * N), "
        for jx in 2:m % N + 1
            body *= "$(rg[jx]), "
        end
        body *= "$(rg[1]), $scr, $zero\n"
        for jx in 3:m % N + 1
            body *= "\tmov\t$(rg[jx]), $((m ÷ N) * N + jx - 1)*8($res)\n"
        end
        body *= "\tmov\t$(rg[1]), $(m + 1)*8($res)\n"
    end

    return pre * body * post * "\n\tret\n"
end

###############################################################################
# 3 ≤ m ≤ 8, 2 ≤ n ≤ m, or m > 8, 2 ≤ n ≤ 8
###############################################################################

function function_body_M(m::Int, n::Int = m)
    if !((3 ≤ m ≤ 8 && 2 ≤ n ≤ m) || (m > 8 && 2 ≤ n ≤ 8))
        error()
    end

    if m ≤ 8
        res_reg = "%rdi"
        ap_reg = "%rsi"
        bp_reg_old = "%rdx"
        bp_reg = "%rcx"
        mov = "\tmov\t$bp_reg_old, $bp_reg\n"
    else
        res_reg = "%rdi"
        ap_reg_old = "%rdx" # is actually bp
        ap_reg = "%rcx"
        bp_reg = "%rsi"
        mov = "\tmov\t$ap_reg_old, $ap_reg\n"
        m, n = n, m
    end

    regs = ["%rax", "%r8", "%r9", "%r10", "%r11", "%rbx", "%rbp", "%r12", "%r13", "%r14", "%r15"]
    preserve_index = 6 # Starting with %rbx, we have to push it to the stack

    numregs = m + 3 # Number of registers used

    regs_perm = deepcopy(regs[1:m + 1])
    # We want the most significant limb to fall on rax
    if n == 2
        regs_perm[1], regs_perm[m + 1] = regs_perm[m + 1], regs_perm[1]
    else
        regs_perm[1], regs_perm[((n - 3) % (m + 1)) + 1] = regs_perm[((n - 3) % (m + 1)) + 1], regs_perm[1]
    end

    scr_reg = regs[numregs - 1]
    zero_reg = regs[numregs]
    zero_reg_32 = reg_32_bit(zero_reg)

    pre = mov
    pre *= "\tmov\t0*8($bp_reg), %rdx\n"
    body = ""
    post = ""

    for jx in preserve_index:numregs
        pre = pre * "\tpush\t$(regs[jx])\n"
        post = "\tpop\t$(regs[jx])\n" * post
    end

    body *= "\n"
    body *= "\txor\t$zero_reg_32, $zero_reg_32\n"
    body *= "\n"

    body *= "\tm$m\t$res_reg, 0, $ap_reg, 0, "
    for jx in 1:m
        body *= "$(regs_perm[jx]), "
    end
    body *= "$(regs_perm[m + 1]), $scr_reg, $zero_reg\n\n"

    for jx in 1:n - 1
        body *= "\tmov\t$jx*8($bp_reg), %rdx\n"
        body *= "\tam$m\t$res_reg, $jx, $ap_reg, 0, "
        for kx in 1:m + 1
            body *= "$(regs_perm[kx]), "
        end
        body *= "$scr_reg, $zero_reg\n"

        # Reorder registers
        regs_perm[1:m], regs_perm[m + 1] = regs_perm[2:m + 1], regs_perm[1]
    end

    body *= "\n"
    for jx in 1:m
        body *= "\tmov\t$(regs_perm[jx]), $(n + jx - 1)*8($res_reg)\n"
    end
    body *= "\n"

    return pre * body * post * "\n\tret\n"
end

###############################################################################
# sqr
###############################################################################

# NOTE: Although this could be generally programmed, just like the mul case,
# we can skip some instructions and additional storage when hardcoding each
# case.

function function_body_sqr_1()
    regs = ["%rcx", "%rax"]

    r0 = regs[1]
    r1 = regs[2] # Important that r1 is rax

    res_reg = "%rdi"
    ap_reg = "%rsi"

    body = ""

    body *= "\tmov\t0*8($ap_reg), %rdx\n"
    body *= "\tmulx\t%rdx, $r0, $r1\n"
    body *= "\tmov\t$r0, 0*8($res_reg)\n"
    body *= "\tmov\t$r1, 1*8($res_reg)\n"

    return body * "\n\tret\n"
end

function function_body_sqr_2()
    n = 2

    res = _regs[1]
    ap = _regs[2]

    this_regs = [_regs[4:end - 1]; __regs]

    w = this_regs[1:3]

    s1 = this_regs[4]
    s2 = _regs[end]

    pre = "\tmov\t0*8($ap), %rdx\n"
    post = ""

    body = ""

    # Calculate upper triangle
    body *= "\txor\t$(reg_32_bit(w[3])), $(reg_32_bit(w[3]))\n"
    body *= "\tmulx\t1*8($ap), $(w[1]), $(w[2])\n"
    body *= "\tadd\t$(w[1]), $(w[1])\n"
    body *= "\tadc\t$(w[2]), $(w[2])\n"
    body *= "\tadc\t$(w[3]), $(w[3])\n"

    # Calculate diagonal and put into res
    body *= "\tmulx\t%rdx, $s1, $s2\n"
    body *= "\tmov\t$s1, 0*8($res)\n"
    body *= "\tadd\t$s2, $(w[1])\n"
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmov\t$(w[1]), 1*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n"
    body *= "\tadc\t$s1, $(w[2])\n"
    body *= "\tadc\t$(w[3]), $s2\n"
    body *= "\tmov\t$(w[2]), 2*8($res)\n"
    body *= "\tmov\t$s2, 3*8($res)\n"

    return pre * body * post * "\n\tret\n"
end

function function_body_sqr_3()
    n = 3

    res = _regs[1]
    ap = _regs[2]

    this_regs = [_regs[4:end - 1]; __regs]

    w = this_regs[1:5]

    scr = _regs[end]

    pre = "\tmov\t0*8($ap), %rdx\n"
    post = ""

    body = ""

    # Calculate upper triangle
    body *= "\txor\t$(reg_32_bit(w[5])), $(reg_32_bit(w[5]))\n"
    body *= "\tmulx\t1*8($ap), $(w[1]), $scr\n" # a0 a1
    body *= "\tmulx\t2*8($ap), $(w[2]), $(w[3])\n" # a0 a2
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tadd\t$scr, $(w[2])\n"
    body *= "\tmulx\t2*8($ap), $scr, $(w[4])\n" # a1 a2
    body *= "\tadc\t$scr, $(w[3])\n"
    body *= "\tadc\t$(w[5]), $(w[4])\n"

    # Double upper triangle
    body *= "\tmov\t0*8($ap), %rdx\n"
    body *= "\tadd\t$(w[1]), $(w[1])\n"
    body *= "\tadc\t$(w[2]), $(w[2])\n"
    body *= "\tadc\t$(w[3]), $(w[3])\n"
    body *= "\tadc\t$(w[4]), $(w[4])\n"
    body *= "\tadc\t$(w[5]), $(w[5])\n"

    # Calculate diagonal and put into res
    body *= "\tmulx\t%rdx, %rdx, $scr\n" # a0^2
    body *= "\tmov\t%rdx, 0*8($res)\n"
    body *= "\tadd\t$scr, $(w[1])\n"
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmov\t$(w[1]), 1*8($res)\n"
    body *= "\tmulx\t%rdx, %rdx, $scr\n" # a1^2
    body *= "\tadc\t%rdx, $(w[2])\n"
    body *= "\tadc\t$scr, $(w[3])\n"
    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tmov\t$(w[2]), 2*8($res)\n"
    body *= "\tmov\t$(w[3]), 3*8($res)\n"
    body *= "\tmulx\t%rdx, %rdx, $scr\n" # a2^2
    body *= "\tadc\t%rdx, $(w[4])\n"
    body *= "\tadc\t$(w[5]), $scr\n"
    body *= "\tmov\t$(w[4]), 4*8($res)\n"
    body *= "\tmov\t$scr, 5*8($res)\n"

    return pre * body * post * "\n\tret\n"
end

function function_body_sqr_4()
    n = 4

    res = _regs[1]
    ap = _regs[2]

    this_regs = [_regs[4:end - 1]; __regs]

    w = this_regs[1:2 * n - 1]

    s1 = "error s1"
    s2 = _regs[end] # Important that this is rax

    pre = "\tmov\t0*8($ap), %rdx\n"
    post = ""
    for jx in 1:2 * n - 1 - length(_regs[4:end - 1])
        pre = pre * "\tpush\t$(__regs[jx])\n"
        post = "\tpop\t$(__regs[jx])\n" * post
    end

    body = ""

    # Calculate upper triangle
    body *= "\txor\t$(reg_32_bit(w[7])), $(reg_32_bit(w[7]))\n"
    body *= "\tmulx\t1*8($ap), $(w[1]), $(w[5])\n" # a0 a1
    body *= "\tmulx\t2*8($ap), $(w[2]), $(w[6])\n" # a0 a2
    body *= "\tmulx\t3*8($ap), $(w[3]), $(w[4])\n" # a0 a3
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tadox\t$(w[5]), $(w[2])\n" # ADOX
    body *= "\tadox\t$(w[6]), $(w[3])\n"
    body *= "\tmulx\t2*8($ap), $(w[5]), $(w[6])\n" # a1 a2
    body *= "\tadcx\t$(w[5]), $(w[3])\n" # ADCX
    body *= "\tadcx\t$(w[6]), $(w[4])\n"
    body *= "\tmulx\t3*8($ap), $(w[6]), $(w[5])\n" # a1 a3
    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tadox\t$(w[6]), $(w[4])\n" # ADOX
    body *= "\tadox\t$(w[7]), $(w[5])\n" # Add zero
    body *= "\tmulx\t3*8($ap), %rdx, $(w[6])\n" # a2 a3
    body *= "\tadcx\t%rdx, $(w[5])\n" # ADCX
    body *= "\tadc\t$(w[7]), $(w[6])\n" # Add zero

    # Double upper triangle
    body *= "\tmov\t0*8($ap), %rdx\n"
    body *= "\tadd\t$(w[1]), $(w[1])\n"
    body *= "\tadc\t$(w[2]), $(w[2])\n"
    body *= "\tadc\t$(w[3]), $(w[3])\n"
    body *= "\tadc\t$(w[4]), $(w[4])\n"
    body *= "\tadc\t$(w[5]), $(w[5])\n"
    body *= "\tadc\t$(w[6]), $(w[6])\n"
    body *= "\tsetc\t$(reg_8_bit(w[7]))\n"

    # Calculate diagonal and put into res
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a0^2
    body *= "\tmov\t%rdx, 0*8($res)\n"
    body *= "\tadd\t$s2, $(w[1])\n"
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmov\t$(w[1]), 1*8($res)\n"
    s1 = w[1]
    w[1] = "error w[1]"
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a1^2
    body *= "\tadc\t%rdx, $(w[2])\n"
    body *= "\tadc\t$s2, $(w[3])\n"
    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tmov\t$(w[2]), 2*8($res)\n"
    body *= "\tmov\t$(w[3]), 3*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a2^2
    body *= "\tadc\t$s1, $(w[4])\n"
    body *= "\tadc\t$s2, $(w[5])\n"
    body *= "\tmov\t3*8($ap), %rdx\n"
    body *= "\tmov\t$(w[4]), 4*8($res)\n"
    body *= "\tmov\t$(w[5]), 5*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a3^2
    body *= "\tadc\t$s1, $(w[6])\n"
    body *= "\tadc\t$(w[7]), $s2\n"
    body *= "\tmov\t$(w[6]), 6*8($res)\n"
    body *= "\tmov\t$s2, 7*8($res)\n"

    return pre * body * post * "\n\tret\n"
end

function function_body_sqr_5()
    n = 5

    res = _regs[1]
    ap = _regs[2]

    this_regs = [_regs[4:end - 1]; __regs]

    w = this_regs[1:2 * n - 1]

    s1 = ""
    s2 = _regs[end] # Important that this is rax

    pre = "\tmov\t0*8($ap), %rdx\n"
    post = ""
    for jx in 1:2 * n - 1 - length(_regs[4:end - 1])
        pre = pre * "\tpush\t$(__regs[jx])\n"
        post = "\tpop\t$(__regs[jx])\n" * post
    end

    body = ""

    # Calculate upper triangle
    body *= "\txor\t$(reg_32_bit(w[9])), $(reg_32_bit(w[9]))\n"
    body *= "\tmulx\t1*8($ap), $(w[1]), $(w[6])\n" # a0 a1
    body *= "\tmulx\t2*8($ap), $(w[2]), $(w[7])\n" # a0 a2
    body *= "\tmulx\t3*8($ap), $(w[3]), $(w[8])\n" # a0 a3
    body *= "\tmulx\t4*8($ap), $(w[4]), $(w[5])\n" # a0 a4
    body *= "\tadcx\t$(w[6]), $(w[2])\n"
    body *= "\tadcx\t$(w[7]), $(w[3])\n"
    body *= "\tadcx\t$(w[8]), $(w[4])\n"
    body *= "\tadcx\t$(w[9]), $(w[5])\n" # Add zero

    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmulx\t2*8($ap), $(w[8]), $(w[6])\n" # a1 a2
    body *= "\tadcx\t$(w[8]), $(w[3])\n" # ADCX
    body *= "\tadcx\t$(w[6]), $(w[4])\n"
    body *= "\tmulx\t3*8($ap), $(w[8]), $(w[6])\n" # a1 a3
    body *= "\tadox\t$(w[8]), $(w[4])\n" # ADOX
    body *= "\tadox\t$(w[6]), $(w[5])\n"
    body *= "\tmulx\t4*8($ap), $(w[8]), $(w[6])\n" # a1 a4
    body *= "\tadcx\t$(w[8]), $(w[5])\n" # ADCX
    body *= "\tadcx\t$(w[9]), $(w[6])\n" # Add zero

    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tmulx\t3*8($ap), $(w[8]), $(w[7])\n" # a2 a3
    body *= "\tadcx\t$(w[8]), $(w[5])\n" # ADCX
    body *= "\tadcx\t$(w[7]), $(w[6])\n"
    body *= "\tmulx\t4*8($ap), $(w[8]), $(w[7])\n" # a2 a4
    body *= "\tadox\t$(w[8]), $(w[6])\n" # ADOX
    body *= "\tadox\t$(w[9]), $(w[7])\n" # Add zero

    body *= "\tmov\t3*8($ap), %rdx\n"
    body *= "\tmulx\t4*8($ap), %rdx, $(w[8])\n" # a3 a4
    body *= "\tadcx\t%rdx, $(w[7])\n" # ADCX
    body *= "\tadc\t$(w[9]), $(w[8])\n" # Add zero

    # Double upper triangle
    body *= "\tmov\t0*8($ap), %rdx\n"
    body *= "\tadd\t$(w[1]), $(w[1])\n"
    body *= "\tadc\t$(w[2]), $(w[2])\n"
    body *= "\tadc\t$(w[3]), $(w[3])\n"
    body *= "\tadc\t$(w[4]), $(w[4])\n"
    body *= "\tadc\t$(w[5]), $(w[5])\n"
    body *= "\tadc\t$(w[6]), $(w[6])\n"
    body *= "\tadc\t$(w[7]), $(w[7])\n"
    body *= "\tadc\t$(w[8]), $(w[8])\n"
    body *= "\tadc\t$(w[9]), $(w[9])\n"

    # Calculate diagonal and put into res
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a0^2
    body *= "\tmov\t%rdx, 0*8($res)\n"
    body *= "\tadd\t$s2, $(w[1])\n"
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmov\t$(w[1]), 1*8($res)\n"
    s1 = w[1]
    w[1] = "error w[1]"
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a1^2
    body *= "\tadc\t%rdx, $(w[2])\n"
    body *= "\tadc\t$s2, $(w[3])\n"
    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tmov\t$(w[2]), 2*8($res)\n"
    body *= "\tmov\t$(w[3]), 3*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a2^2
    body *= "\tadc\t$s1, $(w[4])\n"
    body *= "\tadc\t$s2, $(w[5])\n"
    body *= "\tmov\t3*8($ap), %rdx\n"
    body *= "\tmov\t$(w[4]), 4*8($res)\n"
    body *= "\tmov\t$(w[5]), 5*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a3^2
    body *= "\tadc\t$s1, $(w[6])\n"
    body *= "\tadc\t$s2, $(w[7])\n"
    body *= "\tmov\t4*8($ap), %rdx\n"
    body *= "\tmov\t$(w[6]), 6*8($res)\n"
    body *= "\tmov\t$(w[7]), 7*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a4^2
    body *= "\tadc\t$s1, $(w[8])\n"
    body *= "\tadc\t$(w[9]), $s2\n"
    body *= "\tmov\t$(w[8]), 8*8($res)\n"
    body *= "\tmov\t$s2, 9*8($res)\n"

    return pre * body * post * "\n\tret\n"
end

function function_body_sqr_6()
    n = 6

    res = _regs[1]
    ap = _regs[2]

    this_regs = [_regs[4:end - 1]; __regs]

    w = this_regs[1:2 * n - 1]

    s1 = "error s1"
    s2 = _regs[end] # Important that this is rax

    pre = "\tmov\t0*8($ap), %rdx\n"
    post = ""
    for jx in 1:6
        pre = pre * "\tpush\t$(__regs[jx])\n"
        post = "\tpop\t$(__regs[jx])\n" * post
    end

    body = ""

    # Calculate upper triangle
    body *= "\txor\t$(reg_32_bit(w[11])), $(reg_32_bit(w[11]))\n"
    body *= "\tmulx\t1*8($ap), $(w[1]), $(w[7])\n" # a0 a1
    body *= "\tmulx\t2*8($ap), $(w[2]), $(w[8])\n" # a0 a2
    body *= "\tmulx\t3*8($ap), $(w[3]), $(w[9])\n" # a0 a3
    body *= "\tmulx\t4*8($ap), $(w[4]), $(w[10])\n" # a0 a4
    body *= "\tmulx\t5*8($ap), $(w[5]), $(w[6])\n" # a0 a5
    body *= "\tadcx\t$(w[7]), $(w[2])\n"
    body *= "\tadcx\t$(w[8]), $(w[3])\n"
    body *= "\tadcx\t$(w[9]), $(w[4])\n"
    body *= "\tadcx\t$(w[10]), $(w[5])\n"
    body *= "\tadcx\t$(w[11]), $(w[6])\n" # Add zero

    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmulx\t2*8($ap), $(w[10]), $(w[7])\n" # a1 a2
    body *= "\tadcx\t$(w[10]), $(w[3])\n" # ADCX
    body *= "\tadcx\t$(w[7]), $(w[4])\n"
    body *= "\tmulx\t3*8($ap), $(w[10]), $(w[7])\n" # a1 a3
    body *= "\tadox\t$(w[10]), $(w[4])\n" # ADOX
    body *= "\tadox\t$(w[7]), $(w[5])\n"
    body *= "\tmulx\t4*8($ap), $(w[10]), $(w[7])\n" # a1 a4
    body *= "\tadcx\t$(w[10]), $(w[5])\n" # ADCX
    body *= "\tadcx\t$(w[7]), $(w[6])\n"
    body *= "\tmulx\t5*8($ap), $(w[10]), $(w[7])\n" # a1 a5
    body *= "\tadox\t$(w[10]), $(w[6])\n" # ADOX
    body *= "\tadox\t$(w[11]), $(w[7])\n" # Add zero
    body *= "\tadcx\t$(w[11]), $(w[7])\n" # Add zero

    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tmulx\t3*8($ap), $(w[10]), $(w[8])\n" # a2 a3
    body *= "\tadcx\t$(w[10]), $(w[5])\n" # ADCX
    body *= "\tadcx\t$(w[8]), $(w[6])\n"
    body *= "\tmulx\t4*8($ap), $(w[10]), $(w[8])\n" # a2 a4
    body *= "\tadox\t$(w[10]), $(w[6])\n" # ADOX
    body *= "\tadox\t$(w[8]), $(w[7])\n"
    body *= "\tmulx\t5*8($ap), $(w[10]), $(w[8])\n" # a2 a5
    body *= "\tadcx\t$(w[10]), $(w[7])\n" # ADCX
    body *= "\tadcx\t$(w[11]), $(w[8])\n" # Add zero

    body *= "\tmov\t3*8($ap), %rdx\n"
    body *= "\tmulx\t4*8($ap), $(w[10]), $(w[9])\n" # a3 a4
    body *= "\tadcx\t$(w[10]), $(w[7])\n" # ADCX
    body *= "\tadcx\t$(w[9]), $(w[8])\n"
    body *= "\tmulx\t5*8($ap), $(w[10]), $(w[9])\n" # a3 a5
    body *= "\tadox\t$(w[10]), $(w[8])\n" # ADOX
    body *= "\tadox\t$(w[11]), $(w[9])\n" # Add zero

    body *= "\tmov\t4*8($ap), %rdx\n"
    body *= "\tmulx\t5*8($ap), %rdx, $(w[10])\n" # a4 a5
    body *= "\tadcx\t%rdx, $(w[9])\n" # ADCX
    body *= "\tadc\t$(w[11]), $(w[10])\n" # Add zero

    # Double upper triangle
    body *= "\tmov\t0*8($ap), %rdx\n"
    body *= "\tadd\t$(w[1]), $(w[1])\n"
    body *= "\tadc\t$(w[2]), $(w[2])\n"
    body *= "\tadc\t$(w[3]), $(w[3])\n"
    body *= "\tadc\t$(w[4]), $(w[4])\n"
    body *= "\tadc\t$(w[5]), $(w[5])\n"
    body *= "\tadc\t$(w[6]), $(w[6])\n"
    body *= "\tadc\t$(w[7]), $(w[7])\n"
    body *= "\tadc\t$(w[8]), $(w[8])\n"
    body *= "\tadc\t$(w[9]), $(w[9])\n"
    body *= "\tadc\t$(w[10]), $(w[10])\n"
    body *= "\tsetc\t$(reg_8_bit(w[11]))\n"

    # Calculate diagonal and put into res
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a0^2
    body *= "\tmov\t%rdx, 0*8($res)\n"
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tadd\t$s2, $(w[1])\n"
    body *= "\tmov\t$(w[1]), 1*8($res)\n"
    s1 = w[1]
    w[1] = "error w[1]"
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a1^2
    body *= "\tadc\t%rdx, $(w[2])\n"
    body *= "\tadc\t$s2, $(w[3])\n"
    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tmov\t$(w[2]), 2*8($res)\n"
    body *= "\tmov\t$(w[3]), 3*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a2^2
    body *= "\tadc\t$s1, $(w[4])\n"
    body *= "\tadc\t$s2, $(w[5])\n"
    body *= "\tmov\t3*8($ap), %rdx\n"
    body *= "\tmov\t$(w[4]), 4*8($res)\n"
    body *= "\tmov\t$(w[5]), 5*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a3^2
    body *= "\tadc\t$s1, $(w[6])\n"
    body *= "\tadc\t$s2, $(w[7])\n"
    body *= "\tmov\t4*8($ap), %rdx\n"
    body *= "\tmov\t$(w[6]), 6*8($res)\n"
    body *= "\tmov\t$(w[7]), 7*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a4^2
    body *= "\tadc\t$s1, $(w[8])\n"
    body *= "\tadc\t$s2, $(w[9])\n"
    body *= "\tmov\t5*8($ap), %rdx\n"
    body *= "\tmov\t$(w[8]), 8*8($res)\n"
    body *= "\tmov\t$(w[9]), 9*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a5^2
    body *= "\tadc\t$s1, $(w[10])\n"
    body *= "\tadc\t$(w[11]), $s2\n"
    body *= "\tmov\t$(w[10]), 10*8($res)\n"
    body *= "\tmov\t$s2, 11*8($res)\n"

    return pre * body * post * "\n\tret\n"
end

function function_body_sqr_7()
    n = 7

    res = _regs[1]
    ap = _regs[2]

    this_regs = [_regs[4:end - 1]; __regs]

    w = [this_regs[1:3]; _regs[end]; this_regs[5:11]; ["error w[$ix]" for ix in 12:2 * n - 2]]

    s1 = "error s1"
    s2 = "error s2"
    zero = this_regs[4]

    pre = "\tmov\t0*8($ap), %rdx\n"
    post = ""
    for jx in 1:6
        pre = pre * "\tpush\t$(__regs[jx])\n"
        post = "\tpop\t$(__regs[jx])\n" * post
    end

    body = ""

    # Calculate upper triangle
    body *= "\txor\t$(reg_32_bit(zero)), $(reg_32_bit(zero))\n"
    body *= "\tmulx\t1*8($ap), $(w[1]), $(w[5])\n" # a0 a1
    body *= "\tmulx\t2*8($ap), $(w[2]), $(w[6])\n" # a0 a2
    body *= "\tmulx\t3*8($ap), $(w[3]), $(w[10])\n" # a0 a3
    body *= "\tadcx\t$(w[5]), $(w[2])\n"
    body *= "\tadcx\t$(w[6]), $(w[3])\n"
    body *= "\tmulx\t4*8($ap), $(w[4]), $(w[8])\n" # a0 a4
    body *= "\tmulx\t5*8($ap), $(w[5]), $(w[9])\n" # a0 a5
    body *= "\tmulx\t6*8($ap), $(w[6]), $(w[7])\n" # a0 a6
    body *= "\tadcx\t$(w[10]), $(w[4])\n"
    body *= "\tadcx\t$(w[8]), $(w[5])\n"
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tadcx\t$(w[9]), $(w[6])\n"
    body *= "\tadcx\t$zero, $(w[7])\n"

    body *= "\tmulx\t2*8($ap), $(w[10]), $(w[8])\n" # a1 a2
    body *= "\tmulx\t3*8($ap), $(w[11]), $(w[9])\n" # a1 a3
    body *= "\tadcx\t$(w[10]), $(w[3])\n" # ADCX
    body *= "\tadcx\t$(w[8]), $(w[4])\n"
    body *= "\tadox\t$(w[11]), $(w[4])\n" # ADOX
    body *= "\tadox\t$(w[9]), $(w[5])\n"
    body *= "\tmulx\t4*8($ap), $(w[10]), $(w[8])\n" # a1 a4
    body *= "\tmulx\t5*8($ap), $(w[11]), $(w[9])\n" # a1 a5
    body *= "\tadcx\t$(w[10]), $(w[5])\n" # ADCX
    body *= "\tadcx\t$(w[8]), $(w[6])\n"
    body *= "\tmulx\t6*8($ap), $(w[10]), $(w[8])\n" # a1 a6
    body *= "\tadox\t$(w[11]), $(w[6])\n" # ADOX
    body *= "\tadox\t$(w[9]), $(w[7])\n"
    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tadcx\t$(w[10]), $(w[7])\n" # ADCX
    body *= "\tadox\t$zero, $(w[8])\n"
    body *= "\tadcx\t$zero, $(w[8])\n"

    body *= "\tmulx\t3*8($ap), $(w[11]), $(w[9])\n" # a2 a3
    body *= "\tadox\t$(w[11]), $(w[5])\n" # ADOX
    body *= "\tadox\t$(w[9]), $(w[6])\n"
    body *= "\tmulx\t4*8($ap), $(w[11]), $(w[9])\n" # a2 a4
    body *= "\tadcx\t$(w[11]), $(w[6])\n" # ADCX
    body *= "\tadcx\t$(w[9]), $(w[7])\n"
    body *= "\tmulx\t5*8($ap), $(w[11]), $(w[9])\n" # a2 a5
    body *= "\tadox\t$(w[11]), $(w[7])\n" # ADOX
    body *= "\tadox\t$(w[9]), $(w[8])\n"
    body *= "\tmulx\t6*8($ap), $(w[11]), $(w[9])\n" # a2 a6
    body *= "\tadcx\t$(w[11]), $(w[8])\n" # ADCX
    body *= "\tadox\t$zero, $(w[9])\n"
    body *= "\tmov\t3*8($ap), %rdx\n"
    body *= "\tadcx\t$zero, $(w[9])\n"

    body *= "\tmulx\t4*8($ap), $(w[11]), $(w[10])\n" # a3 a4
    body *= "\tadcx\t$(w[11]), $(w[7])\n" # ADCX
    body *= "\tadcx\t$(w[10]), $(w[8])\n"
    body *= "\tmulx\t5*8($ap), $(w[11]), $(w[10])\n" # a3 a5
    body *= "\tadox\t$(w[11]), $(w[8])\n" # ADOX
    body *= "\tadox\t$(w[10]), $(w[9])\n"
    body *= "\tmulx\t6*8($ap), $(w[11]), $(w[10])\n" # a3 a6
    body *= "\tmov\t4*8($ap), %rdx\n"
    body *= "\tadcx\t$(w[11]), $(w[9])\n" # ADCX
    body *= "\tadcx\t$zero, $(w[10])\n"

    w[12] = zero
    zero = "error zero"
    body *= "\tmulx\t5*8($ap), $(w[12]), $(w[11])\n" # a4 a5
    body *= "\tadcx\t$(w[12]), $(w[9])\n" # ADCX
    body *= "\tadcx\t$(w[11]), $(w[10])\n"
    body *= "\tmulx\t6*8($ap), $(w[12]), $(w[11])\n" # a4 a6
    body *= "\tmov\t\$0, %edx\n"
    body *= "\tadox\t$(w[12]), $(w[10])\n" # ADOX
    body *= "\tadox\t%rdx, $(w[11])\n"
    body *= "\tmov\t5*8($ap), %rdx\n"

    body *= "\tmulx\t6*8($ap), %rdx, $(w[12])\n" # a5 a6
    body *= "\tadcx\t%rdx, $(w[11])\n" # ADCX
    body *= "\tadc\t\$0, $(w[12])\n"
    body *= "\ttest\t%al, %al\n"

    # Double upper triangle
    body *= "\tmov\t0*8($ap), %rdx\n"
    body *= "\tadcx\t$(w[1]), $(w[1])\n"
    body *= "\tadcx\t$(w[2]), $(w[2])\n"
    body *= "\tadcx\t$(w[3]), $(w[3])\n"
    body *= "\tadcx\t$(w[4]), $(w[4])\n"
    body *= "\tadcx\t$(w[5]), $(w[5])\n"
    body *= "\tadcx\t$(w[6]), $(w[6])\n"
    body *= "\tpush\t$(w[4])\n"
    body *= "\tadcx\t$(w[7]), $(w[7])\n"
    body *= "\tadcx\t$(w[8]), $(w[8])\n"
    body *= "\tadcx\t$(w[9]), $(w[9])\n"
    body *= "\tadcx\t$(w[10]), $(w[10])\n"
    body *= "\tadcx\t$(w[11]), $(w[11])\n"
    body *= "\tadcx\t$(w[12]), $(w[12])\n"

    # Calculate diagonal and put into res
    s2 = w[4]
    w[4] = "error w[4]"
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a0^2
    body *= "\tmov\t%rdx, 0*8($res)\n"
    body *= "\tadox\t$s2, $(w[1])\n"
    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmov\t$(w[1]), 1*8($res)\n"
    w[4] = w[1]
    w[1] = "error w[1]"
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a1^2
    body *= "\tadox\t%rdx, $(w[2])\n"
    body *= "\tadox\t$s2, $(w[3])\n"
    body *= "\tpop\t$(w[4])\n"
    body *= "\tmov\t2*8($ap), %rdx\n"
    body *= "\tmov\t$(w[2]), 2*8($res)\n"
    body *= "\tmov\t$(w[3]), 3*8($res)\n"
    s1 = w[2]
    w[2] = "error w[2]"
    body *= "\tmulx\t%rdx, %rdx, $s2\n" # a2^2
    body *= "\tadox\t%rdx, $(w[4])\n"
    body *= "\tadox\t$s2, $(w[5])\n"
    body *= "\tmov\t3*8($ap), %rdx\n"
    body *= "\tmov\t$(w[4]), 4*8($res)\n"
    body *= "\tmov\t$(w[5]), 5*8($res)\n"
    zero = w[4]
    w[4] = "error w[4]"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a3^2
    body *= "\tadox\t$s1, $(w[6])\n"
    body *= "\tadox\t$s2, $(w[7])\n"
    body *= "\tmov\t4*8($ap), %rdx\n"
    body *= "\tmov\t$(w[6]), 6*8($res)\n"
    body *= "\tmov\t$(w[7]), 7*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a4^2
    body *= "\tadox\t$s1, $(w[8])\n"
    body *= "\tadox\t$s2, $(w[9])\n"
    body *= "\tmov\t5*8($ap), %rdx\n"
    body *= "\tmov\t$(w[8]), 8*8($res)\n"
    body *= "\tmov\t$(w[9]), 9*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a5^2
    body *= "\tadox\t$s1, $(w[10])\n"
    body *= "\tadox\t$s2, $(w[11])\n"
    body *= "\tmov\t\$0, $(reg_32_bit(zero))\n"
    body *= "\tmov\t6*8($ap), %rdx\n"
    body *= "\tmov\t$(w[10]), 10*8($res)\n"
    body *= "\tmov\t$(w[11]), 11*8($res)\n"
    body *= "\tmulx\t%rdx, $s1, $s2\n" # a6^2
    body *= "\tadox\t$s1, $(w[12])\n"
    body *= "\tadox\t$zero, $s2\n"
    body *= "\tadcx\t$zero, $s2\n"
    body *= "\tmov\t$(w[12]), 12*8($res)\n"
    body *= "\tmov\t$s2, 13*8($res)\n"

    return pre * body * post * "\n\tret\n"
end

# STRATEGY: We compute (a1 B + a0)^2 = a1^2 B^2 + a0^2 + 2 a0 a1 B by
# computing 2 a0 a1 B into stack, then adding that onto a1^2 B^2 + a0^2.
function function_body_sqr_8()
    n = 8

    res = _regs[1]
    ap = _regs[2]

    this_regs = [_regs[4:end - 1]; __regs]

    w = [this_regs[1:3]; _regs[end]; this_regs[5:11]; ["error w[$ix]" for ix in 12:2 * n - 2]]

    s = ["error s[$ix]" for ix in 1:n]
    zero = this_regs[4]

    pre = "\tmov\t0*8($ap), %rdx\n"
    post = ""
    for jx in 1:6
        pre = pre * "\tpush\t$(__regs[jx])\n"
        post = "\tpop\t$(__regs[jx])\n" * post
    end

    body = ""

    # Calculate upper triangle
    body *= "\txor\t$(reg_32_bit(zero)), $(reg_32_bit(zero))\n"
    body *= "\tmulx\t4*8($ap), $(w[4]), $(s[1])\n" # a0 a4
    body *= "\tmulx\t5*8($ap), $(w[5]), $(s[2])\n" # a0 a5
    body *= "\tmulx\t6*8($ap), $(w[6]), $(s[3])\n" # a0 a6
    body *= "\tmulx\t7*8($ap), $(w[7]), $(w[8])\n" # a0 a7
    body *= "\tadcx\t$(s[1]), $(w[5])\n"
    body *= "\tadcx\t$(s[2]), $(w[6])\n"
    body *= "\tadcx\t$(s[3]), $(w[7])\n"
    body *= "\tadcx\t$(zero), $(w[8])\n"

    body *= "\tmov\t1*8($ap), %rdx\n"
    body *= "\tmulx\t4*8($ap), $(s[1]), $(s[2])\n" # a1 a4
    body *= "\tmulx\t5*8($ap), $(s[3]), $(s[4])\n" # a1 a5
    body *= "\tadox\t$(s[2]), $(s[3])\n"
    body *= "\tadcx\t$(s[1]), $(w[5])\n"
    body *= "\tadcx\t$(s[3]), $(w[6])\n"
    body *= "\tmulx\t6*8($ap), $(s[1]), $(w[10])\n" # a1 a6
    body *= "\tmulx\t7*8($ap), $(w[11]), $(w[9])\n" # a1 a7
    body *= "\tadcx\t$(s[3]), $(w[7])\n"
    body *= "\tadcx\t$(zero), $(w[8])\n"

    return pre * body * post * "\n\tret\n"
end

###############################################################################
# Generate file
###############################################################################

function gen_mul(m::Int, n::Int = m, nofile::Bool = false)
    (pre, post) = function_pre_post("flint_mpn_mul_$(m)_$n")
    macros = ""
    functionbody = ""
    if m == 1 && n == 1
        functionbody = function_body_1(m, n)
    elseif m == 2 && n == 1
        functionbody = function_body_2_1(m, n)
    elseif m == 2 && n == 2
        functionbody = function_body_2(m, n)
    elseif n == 1
        macros = mul_1_macro(m)
        functionbody = function_body_M_1(m)
    elseif m ≤ 8
        macros = mulM_macro(m) * "\n" * addmulM_macro(m) * "\n"
        functionbody = function_body_M(m, n)
    elseif m > 8 && n != 2
        macros = mulM_macro(n) * "\n" * addmulM_macro(n) * "\n"
        functionbody = function_body_M(m, n)
    elseif n == 2
        N = M_2_getN(m)
        macros = mulM_macro(N) * "\n"
        macros *= addmulM_macro(N) * "\n"
        if m ÷ N >= 2
            macros *= mulM_macro(N, chain = true) * "\n"
        end
        if m % N > 1
            macros *= mulM_macro(m % N, chain = true) * "\n"
            macros *= addmulM_macro(m % N) * "\n"
        end
        functionbody = function_body_M_2(m)
    else
        error("This won't work for m = $m and n = $n")
    end

    str = "$copyright\n$preamble\n$macros$pre$functionbody$post"

    if nofile
        print(str)
    else
        path = String(@__DIR__) * "/../src/mpn_extras/broadwell/mul_$(m)_$n.asm"
        file = open(path, "w")
        write(file, str)
        close(file)
    end
end

function gen_sqr(m::Int, nofile::Bool = false)
    (pre, post) = function_pre_post("flint_mpn_sqr_$m")
    functionbody = ""
    if m == 1
        functionbody = function_body_sqr_1()
    elseif m == 2
        functionbody = function_body_sqr_2()
    elseif m == 3
        functionbody = function_body_sqr_3()
    elseif m == 4
        functionbody = function_body_sqr_4()
    elseif m == 5
        functionbody = function_body_sqr_5()
    elseif m == 6
        functionbody = function_body_sqr_6()
    elseif m == 7
        functionbody = function_body_sqr_7()
    elseif m == 8
        println("Warning: Generating flint_mpn_sqr_8 which has about the same speed as GMP.")
        functionbody = function_body_sqr_8()
    else
        error("This won't work for m = $m")
    end

    str = "$copyright\n$preamble\n$pre$functionbody$post"

    if nofile
        print(str)
    else
        path = String(@__DIR__) * "/../src/mpn_extras/broadwell/sqr_$m.asm"
        file = open(path, "w")
        write(file, str)
        close(file)
    end
end

function gen_all()
    for mx in 1:16
        for nx in 1:min(mx, 8)
            gen_mul(mx, nx)
        end
    end

    for mx in 1:7
        gen_sqr(mx)
    end
end
