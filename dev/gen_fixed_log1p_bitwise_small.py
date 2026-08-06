#!/usr/bin/env python3
"""Regenerate src/fixed/log1p_bitwise_rs_small.inc (run from the
top-level FLINT directory: python3 dev/gen_fixed_log1p_bitwise_small.py
> src/fixed/log1p_bitwise_rs_small.inc): per-size register
implementations of the fixed_log1p_bitwise_rs reduction for
n = 3..7, mirroring the exp small path (dev/
gen_fixed_exp_bitwise_small.py).

Everything lives in named scalars: the deficit d0..d{n-1}, the
product fraction p0..p{n-1} (the units limb of P is implicitly 1:
accepted factors satisfy P + v <= X < 2, so the fraction never
carries out), and the accumulated logarithms a0..a{n-1}.  Each
window c (FLINT_BITS-aligned) runs the single-limb model: the fast
path rejects on h < lt (lt being one funnel shift of P's top limbs)
with a stray-bit check for c >= 1 (the deficit obeys only
D < P 2^-i < 2^(1-i), so d{n-c} may hold one bit; when it does the
factor certainly fits).  Everything else -- window boundaries, model
ties, certain accepts and the extra step at i = r -- goes through a
single masked step: one borrow-capturing subtraction chain doubles
as the exact comparison and, through the resulting all-ones/zero
mask, selects the updated deficit and gates the additions into P and
the accumulator.  The extra step's window index depends on the
runtime r, so it is emitted as a switch on r / FLINT_BITS with the
boundary (b == 0) form handled separately in each case."""

SUB = {2: "sub_ddmmss", 3: "sub_dddmmmsss", 4: "sub_ddddmmmmssss",
       5: "sub_dddddmmmmmsssss", 6: "sub_ddddddmmmmmmssssss",
       7: "sub_dddddddmmmmmmmsssssss",
       8: "sub_ddddddddmmmmmmmmssssssss"}
ADD = {2: "add_ssaaaa", 3: "add_sssaaaaaa", 4: "add_ssssaaaaaaaa",
       5: "add_sssssaaaaaaaaaa", 6: "add_ssssssaaaaaaaaaaaa",
       7: "add_sssssssaaaaaaaaaaaaaa",
       8: "add_ssssssssaaaaaaaaaaaaaaaa"}

NMIN, NMAX = 3, 7


def masked_step(o, n, c, vt):
    """one masked accept of the factor with fraction limbs vt[j]
    (j = 0..n-1-c) and, for boundary steps, an implicit leading 1 at
    limb n - c; vt[n - c] is "UWORD(1)" then, else absent.  Emits the
    borrow-chain compare/select and the gated updates; assumes the
    table pointer Lp for the factor is set."""
    top = n - c if (c >= 1 or len(vt) > n - c) else n - 1 - c
    span = top + 1                       # limbs 0..top of D
    dl = ["d%d" % j for j in range(top, -1, -1)]
    vl = []
    for j in range(top, -1, -1):
        vl.append(vt[j] if j < len(vt) else "UWORD(0)")
    el = ["e%d" % j for j in range(top, -1, -1)]
    # borrow-capturing chain: one extra (0 - 0 - borrow) limb yields
    # bw = all-ones exactly when v > D
    o("        %s(bw, %s," % (SUB[span + 1], ", ".join(el)))
    o("            UWORD(0), %s," % ", ".join(dl))
    o("            UWORD(0), %s);" % ", ".join(vl))
    o("        m = ~bw;                /* accept iff no borrow */")
    for j in range(top + 1):
        o("        d%d = (d%d & bw) | (e%d & m);" % (j, j, j))
    # P += v & m: v occupies limbs 0..top, carries may ripple to the
    # fraction top but never out (P + v <= X < 2)
    pl = ["p%d" % j for j in range(n - 1, -1, -1)]
    va = ["(%s) & m" % vt[j] if j < len(vt) else "UWORD(0)"
          for j in range(n - 1, -1, -1)]
    o("        %s(%s," % (ADD[n], ", ".join(pl)))
    o("            %s," % ", ".join(pl))
    o("            %s);" % ", ".join(va))
    # acc += L & m over all n limbs (no carry out: acc < log 2)
    al = ["a%d" % j for j in range(n - 1, -1, -1)]
    la = ["Lp[%d] & m" % j for j in range(n - 1, -1, -1)]
    o("        %s(%s," % (ADD[n], ", ".join(al)))
    o("            %s," % ", ".join(al))
    o("            %s);" % ", ".join(la))


OPT_R = {3: 10, 4: 26, 5: 31, 6: 30}


def gen_one(o, n, opt=False, fn_name=None, series_name=None,
        r_const=None, static=True):
    if fn_name is None:
        fn_name = ("_fixed_log1p_bitwise_rs_opt_%d" % n) if opt \
            else ("_fixed_log1p_bitwise_rs_%d" % n)
    if series_name is None and opt:
        series_name = "_fixed_atanh_rs_opt_%d" % n
    if r_const is None and opt:
        r_const = OPT_R[n]
    o("static void" if static else "void")
    if opt:
        o("%s(nn_ptr res, nn_srcptr x)" % fn_name)
    else:
        o("%s(nn_ptr res, nn_srcptr x, int r)" % fn_name)
    o("{")
    if opt:
        o("    const int r = %d;" % r_const)
    o("    ulong " + ", ".join("d%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("p%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("a%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("e%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("v%d" % j for j in range(n)) + ";")
    o("    ulong bw, m, h, lt;")
    o("    ulong S[%d], nd[%d], t[%d], y[%d];" % (n + 1, 2 * n, n, n))
    o("    slong i, nc;")
    o("")
    lo = 16 if n <= 4 else 32
    if not opt:
        o("    r = FLINT_MAX(r, %d);" % lo)
        o("    r = FLINT_MIN(r, FLINT_BITS * %d - 16);" % n)
    o("")
    o("    _fixed_exp_logs_ensure(%d, r);" % n)
    o("    nc = _fixed_exp_logs_n;")
    o("")
    for j in range(n):
        o("    d%d = x[%d];" % (j, j))
    for j in range(n):
        o("    p%d = 0;" % j)
    for j in range(n):
        o("    a%d = 0;" % j)
    o("")
    o("#define LP(ii) (_fixed_exp_logs + (ii) * nc + (nc - %d))" % n)
    o("")

    cmax = n - 1 if r_const is None else min(n - 1, r_const // 64)
    for c in range(cmax + 1):
        i0 = 1 if c == 0 else 64 * c
        o("    /* window %d */" % c)
        guard = "" if (c == 0 or r_const is not None) \
            else "if (r >= %d)\n    " % (64 * c)
        o("    %s{" % guard)
        o("        nn_srcptr Lp;")
        o("")
        # boundary exact step at i0
        o("        /* boundary step i = %d */" % i0)
        o("        Lp = LP(%d);" % i0)
        if c == 0:
            # i0 = 1: funnel form with b = 1
            vt = []
            for j in range(n - 1):
                vt.append("MPN_RIGHT_SHIFT_LOW(p%d, p%d, 1)"
                          % (j + 1, j))
            vt.append("MPN_RIGHT_SHIFT_LOW(UWORD(1), p%d, 1)" % (n - 1))
            for j in range(n):
                o("        v%d = %s;" % (j, vt[j]))
            masked_step(o, n, 0, ["v%d" % j for j in range(n)])
        else:
            # i0 = 64c: v = (1, p_{n-1}, ..., p_c)
            vt = ["p%d" % (j + c) for j in range(n - c)] + ["UWORD(1)"]
            masked_step(o, n, c, vt)
        o("")
        o("        h = d%d;" % (n - 1 - c))
        o("")
        if r_const is None:
            o("        for (i = %d; i <= FLINT_MIN((slong) r, %d); i++)"
              % (i0 + 1, 64 * c + 63))
        else:
            o("        for (i = %d; i <= %d; i++)"
              % (i0 + 1, min(r_const, 64 * c + 63)))
        o("        {")
        o("            int b = (int) (i - %d);" % (64 * c))
        o("")
        o("            lt = MPN_RIGHT_SHIFT_LOW(UWORD(1), p%d, b);"
          % (n - 1))
        if c >= 1:
            o("            if (h < lt && d%d == 0)" % (n - c))
        else:
            o("            if (h < lt)")
        o("                continue;    /* certain reject */")
        o("")
        o("            Lp = LP(i);")
        for j in range(n - 1 - c):
            o("            v%d = MPN_RIGHT_SHIFT_LOW(p%d, p%d, b);"
              % (j, j + c + 1, j + c))
        o("            v%d = lt;" % (n - 1 - c))
        body = []
        masked_step(lambda ln: body.append(ln), n, c,
                    ["v%d" % j for j in range(n - c)])
        for ln in body:
            o("    " + ln)
        o("            h = d%d;" % (n - 1 - c))
        o("        }")
        o("    }")
        o("")

    # extra step at i = r
    o("    /* one extra step at i = r absorbs the truncation creep */")
    o("    {")
    o("        nn_srcptr Lp = LP((slong) r);")
    if r_const is None:
        o("        int b = (int) (r & (FLINT_BITS - 1));")
        o("")
        o("        switch (r / FLINT_BITS)")
        o("        {")
        for c in range(n):
            o("        case %d:" % c)
            if c >= 1:
                o("            if (b == 0)")
                o("            {")
                vt = ["p%d" % (j + c) for j in range(n - c)] + ["UWORD(1)"]
                body = []
                masked_step(lambda ln: body.append(ln), n, c, vt)
                for ln in body:
                    o("        " + ln)
                o("                break;")
                o("            }")
            for j in range(n - 1 - c):
                o("            v%d = MPN_RIGHT_SHIFT_LOW(p%d, p%d, b);"
                  % (j, j + c + 1, j + c))
            o("            v%d = MPN_RIGHT_SHIFT_LOW(UWORD(1), p%d, b);"
              % (n - 1 - c, n - 1))
            body = []
            masked_step(lambda ln: body.append(ln), n, c,
                        ["v%d" % j for j in range(n - c)])
            for ln in body:
                o("    " + ln)
            o("            break;")
        o("        }")
    else:
        c = r_const // 64
        b = r_const % 64
        if c >= 1 and b == 0:
            vt = ["p%d" % (j + c) for j in range(n - c)] + ["UWORD(1)"]
            masked_step(o, n, c, vt)
        else:
            o("        const int b = %d;" % b)
            o("")
            for j in range(n - 1 - c):
                o("        v%d = MPN_RIGHT_SHIFT_LOW(p%d, p%d, b);"
                  % (j, j + c + 1, j + c))
            o("        v%d = MPN_RIGHT_SHIFT_LOW(UWORD(1), p%d, b);"
              % (n - 1 - c, n - 1))
            masked_step(o, n, c, ["v%d" % j for j in range(n - c)])
    o("    }")
    o("")
    o("#undef LP")
    o("")
    # division
    o("    /* t = D 2^(FLINT_BITS %d) / (X + P) */" % n)
    o("    %s(%s," % (ADD[n + 1],
        ", ".join("S[%d]" % j for j in range(n, -1, -1))))
    o("        UWORD(1), %s," % ", ".join("x[%d]" % j
        for j in range(n - 1, -1, -1)))
    o("        UWORD(1), %s);" % ", ".join("p%d" % j
        for j in range(n - 1, -1, -1)))
    for j in range(n):
        o("    nd[%d] = 0;" % j)
    for j in range(n):
        o("    nd[%d] = d%d;" % (n + j, j))
    o("    flint_mpn_zero(t, %d);" % n)
    o("    if (%s)" % " || ".join("d%d != 0" % j
        for j in range(n - 1, -1, -1)))
    o("        mpn_tdiv_qr(t, nd, 0, nd, %d, S, %d);" % (2 * n, n + 1))
    o("")
    # atanh (r < 32 reaches here only for n <= 4: at n >= 5 the
    # fixed-length 2^-16 series costs more than the reduction saves,
    # so r is clamped to 32 above)
    if opt:
        o("    %s(y, t);" % series_name)
    else:
        o("    fixed_atanh_rs(y, t, %d);" % n)
    o("")
    # res = acc + 2 atanh(t)
    o("    %s(%s," % (ADD[n],
        ", ".join("res[%d]" % j for j in range(n - 1, -1, -1))))
    o("        %s," % ", ".join("a%d" % j for j in range(n - 1, -1, -1)))
    yl = []
    for j in range(n - 1, -1, -1):
        if j == 0:
            yl.append("y[0] << 1")
        else:
            yl.append("(y[%d] << 1) | (y[%d] >> (FLINT_BITS - 1))"
                      % (j, j - 1))
    o("        %s);" % ", ".join(yl))
    o("}")
    o("")



def gen_one_named(n, r, fn_name, series_name, static):
    L = []
    gen_one(L.append, n, opt=True, fn_name=fn_name,
            series_name=series_name, r_const=r, static=static)
    return L

if __name__ == "__main__":
    import sys
    sys.exit("this module is a body-emitter library for "
             "dev/tune_fixed.py; run that instead")
