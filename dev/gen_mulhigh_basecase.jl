#
#   Copyright (C) 2024 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 2.1 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

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

copyright = "#
#   Copyright (C) 2024 Albin Ahlbäck
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
# mulhigh
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
# Generate file
###############################################################################

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

function gen_all()
    for m in 1:12
        gen_mulhigh(m)
    end
end
