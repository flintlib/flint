#!/usr/bin/env python3
"""Regenerate src/fixed/atan_bitwise_rs_small.inc (run from the
top-level FLINT directory: python3 dev/gen_fixed_atan_bitwise_small.py
> src/fixed/atan_bitwise_rs_small.inc): per-size register
implementations of the fixed_atan_bitwise_rs vectoring for n = 1..7,
mirroring the log1p small path (dev/gen_fixed_log1p_bitwise_small.py).

Everything lives in named scalars: the rotating vector's components
X (fraction x0..x{n-1}; the units limb is implicitly 1, since X only
grows, from 1 to at most prod sqrt(1 + 4^-i) < 1.1644, so the
fraction never carries out) and Y (y0..y{n-1}, with 0 <= Y < X), plus
the accumulated angles a0..a{n-1} (their total is below
sum_{i>=1} atan(2^-i) ~ 0.898 < 1, so that sum never carries out
either).

Unlike the log reduction, which applies ONE funnel-shifted vector to
both of its state words, a vectoring step needs two, both read from
the OLD components: v = trunc(X >> i), which is what Y is compared
against and what is subtracted from it, and w = trunc(Y >> i), which
is added to X.  Each window c (FLINT_BITS-aligned) runs the
single-limb model: the fast path rejects on h < lt (lt being one
funnel shift of X's top limbs) with a stray-bit check for c >= 1
(the residual obeys only Y < X 2^-i < 2^(1-i), so y{n-c} may hold one
bit; when it does the factor certainly fits).  Everything else --
window boundaries, model ties, certain accepts and the step at
i = r -- goes through a single masked step: one borrow-capturing
subtraction chain doubles as the exact comparison and, through the
resulting all-ones/zero mask, selects the updated Y and gates the
additions into X and the accumulator.  The step at i = r is emitted
as a switch on r / FLINT_BITS (the boundary b == 0 form handled
separately in each case) and run TWICE: entering it Y < X 2^-(r-1),
so at most two subtractions of trunc(X >> r) are needed to restore
Y < trunc(X >> r) <= X 2^-r, the series contract.
"""

SUB = {2: "sub_ddmmss", 3: "sub_dddmmmsss", 4: "sub_ddddmmmmssss",
       5: "sub_dddddmmmmmsssss", 6: "sub_ddddddmmmmmmssssss",
       7: "sub_dddddddmmmmmmmsssssss",
       8: "sub_ddddddddmmmmmmmmssssssss",
       9: "sub_dddddddddmmmmmmmmmsssssssss"}
ADD = {1: None, 2: "add_ssaaaa", 3: "add_sssaaaaaa",
       4: "add_ssssaaaaaaaa", 5: "add_sssssaaaaaaaaaa",
       6: "add_ssssssaaaaaaaaaaaa",
       7: "add_sssssssaaaaaaaaaaaaaa",
       8: "add_ssssssssaaaaaaaaaaaaaaaa"}

NMIN, NMAX = 1, 7


def add_into(o, ind, dst, src, n):
    """dst[0..n-1] += src[0..n-1], no carry out"""
    if n == 1:
        o("%s%s0 += %s;" % (ind, dst, src[0]))
        return
    dl = ["%s%d" % (dst, j) for j in range(n - 1, -1, -1)]
    sl = [src[j] for j in range(n - 1, -1, -1)]
    o("%s%s(%s," % (ind, ADD[n], ", ".join(dl)))
    o("%s    %s," % (ind, ", ".join(dl)))
    o("%s    %s);" % (ind, ", ".join(sl)))


def masked_step(o, n, c, vt, wt, ind="        "):
    """one masked accept of the factor at index i: vt holds the
    fraction limbs of trunc(X >> i) (with an implicit leading 1 at
    limb n - c for the b == 0 boundary form), wt those of
    trunc(Y >> i).  Emits the borrow-chain compare/select on Y and
    the gated updates of X and the accumulator; assumes Ap points at
    the tabulated angle."""
    top = n - c if (c >= 1 and len(vt) > n - 1 - c) else n - 1 - c
    span = top + 1
    yl = ["y%d" % j for j in range(top, -1, -1)]
    vl = [(vt[j] if j < len(vt) else "UWORD(0)")
          for j in range(top, -1, -1)]
    el = ["e%d" % j for j in range(top, -1, -1)]
    # borrow-capturing chain: one extra (0 - 0 - borrow) limb yields
    # bw = all-ones exactly when v > Y
    o("%s%s(bw, %s," % (ind, SUB[span + 1], ", ".join(el)))
    o("%s    UWORD(0), %s," % (ind, ", ".join(yl)))
    o("%s    UWORD(0), %s);" % (ind, ", ".join(vl)))
    o("%sm = ~bw;                /* accept iff no borrow */" % ind)
    # X and the accumulator are updated BEFORE Y is overwritten: both
    # vectors of a vectoring step read the OLD components, and at a
    # b == 0 boundary wt names the y registers themselves rather than
    # precomputed temporaries.
    # X += w & m (the fraction cannot carry out: X < 1.1644)
    wa = [("(%s) & m" % wt[j]) if j < len(wt) else "UWORD(0)"
          for j in range(n)]
    add_into(o, ind, "x", wa, n)
    # acc += A_i & m (the total stays below 0.898)
    aa = ["Ap[%d] & m" % j for j in range(n)]
    add_into(o, ind, "a", aa, n)
    for j in range(top + 1):
        o("%sy%d = (y%d & bw) | (e%d & m);" % (ind, j, j, j))


def vw_shifted(o, n, c, ind):
    """emit v0.. and w0.. for a runtime shift b within window c"""
    for j in range(n - 1 - c):
        o("%sv%d = MPN_RIGHT_SHIFT_LOW(x%d, x%d, b);"
          % (ind, j, j + c + 1, j + c))
    o("%sv%d = lt;" % (ind, n - 1 - c))
    for j in range(n - 1 - c):
        o("%sw%d = MPN_RIGHT_SHIFT_LOW(y%d, y%d, b);"
          % (ind, j, j + c + 1, j + c))
    o("%sw%d = MPN_RIGHT_SHIFT_LOW(UWORD(0), y%d, b);"
      % (ind, n - 1 - c, n - 1))


# Fully specialized sizes: r is a compile-time constant and the series
# is the hand-written one built for exactly that r (see
# src/fixed/trig_rs_opt_hand.inc).  The hand series is so much cheaper
# than the generated one that the optimal r falls sharply.
OPT_R = {1: 4, 2: 6, 3: 22, 4: 20, 5: 18, 6: 20, 7: 19}


def gen_one(o, n, opt=False, fn_name=None, series_name=None,
        r_const=None, static=True):
    if fn_name is None:
        fn_name = ("_fixed_atan_bitwise_rs_opt_%d" % n) if opt \
            else ("_fixed_atan_bitwise_rs_%d" % n)
    if series_name is None and opt:
        series_name = "_fixed_atan_rs_opt_%d" % n
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
    o("    ulong " + ", ".join("x%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("y%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("a%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("e%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("v%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("w%d" % j for j in range(n)) + ";")
    o("    ulong bw, m, h, lt;")
    o("    ulong S[%d], nd[%d], t[%d], q[%d];" % (n + 1, 2 * n, n, n))
    o("    slong i, nc;")
    o("    int nz;")
    o("")
    if not opt:
        o("    r = FLINT_MAX(r, 16);")
        o("    r = FLINT_MIN(r, FLINT_BITS * %d - 16);" % n)
    o("")
    o("    _fixed_atans_ensure(%d, r);" % n)
    o("    nc = _fixed_atans_n;")
    o("")
    for j in range(n):
        o("    x%d = 0;" % j)
    for j in range(n):
        o("    y%d = x[%d];" % (j, j))
    for j in range(n):
        o("    a%d = 0;" % j)
    o("")
    o("#define AP(ii) (_fixed_atans + (ii) * nc + (nc - %d))" % n)
    o("")

    cmax = n - 1 if r_const is None else min(n - 1, r_const // 64)
    for c in range(cmax + 1):
        i0 = 1 if c == 0 else 64 * c
        o("    /* window %d */" % c)
        guard = "" if (c == 0 or r_const is not None) \
            else "if (r >= %d)\n    " % (64 * c)
        o("    %s{" % guard)
        o("        nn_srcptr Ap;")
        o("")
        o("        /* boundary step i = %d */" % i0)
        o("        Ap = AP(%d);" % i0)
        if c == 0:
            # i0 = 1: funnel form with b = 1, X's units limb is 1
            for j in range(n - 1):
                o("        v%d = MPN_RIGHT_SHIFT_LOW(x%d, x%d, 1);"
                  % (j, j + 1, j))
            o("        v%d = MPN_RIGHT_SHIFT_LOW(UWORD(1), x%d, 1);"
              % (n - 1, n - 1))
            for j in range(n - 1):
                o("        w%d = MPN_RIGHT_SHIFT_LOW(y%d, y%d, 1);"
                  % (j, j + 1, j))
            o("        w%d = MPN_RIGHT_SHIFT_LOW(UWORD(0), y%d, 1);"
              % (n - 1, n - 1))
            masked_step(o, n, 0, ["v%d" % j for j in range(n)],
                        ["w%d" % j for j in range(n)])
        else:
            # i0 = 64c, b = 0: v = (x_c, .., x_{n-1}, 1), w = (y_c, ..)
            vt = ["x%d" % (j + c) for j in range(n - c)] + ["UWORD(1)"]
            wt = ["y%d" % (j + c) for j in range(n - c)]
            masked_step(o, n, c, vt, wt)
        o("")
        o("        h = y%d;" % (n - 1 - c))
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
        o("            lt = MPN_RIGHT_SHIFT_LOW(UWORD(1), x%d, b);"
          % (n - 1))
        if c >= 1:
            o("            if (h < lt && y%d == 0)" % (n - c))
        else:
            o("            if (h < lt)")
        o("                continue;    /* certain reject */")
        o("")
        o("            Ap = AP(i);")
        vw_shifted(o, n, c, "            ")
        masked_step(o, n, c, ["v%d" % j for j in range(n - c)],
                    ["w%d" % j for j in range(n - c)],
                    ind="            ")
        o("            h = y%d;" % (n - 1 - c))
        o("        }")
        o("    }")
        o("")

    # the step at i = r, run twice (see the module comment)
    o("    /* the step at i = r restores Y < trunc(X >> r) <= X 2^-r;")
    o("       entering it Y < X 2^-(r-1), so two passes suffice */")
    o("    for (nz = 0; nz < 2; nz++)")
    o("    {")
    o("        nn_srcptr Ap = AP((slong) r);")
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
                vt = ["x%d" % (j + c) for j in range(n - c)] + ["UWORD(1)"]
                wt = ["y%d" % (j + c) for j in range(n - c)]
                masked_step(o, n, c, vt, wt, ind="                ")
                o("                break;")
                o("            }")
            o("            lt = MPN_RIGHT_SHIFT_LOW(UWORD(1), x%d, b);"
              % (n - 1))
            vw_shifted(o, n, c, "            ")
            masked_step(o, n, c, ["v%d" % j for j in range(n - c)],
                        ["w%d" % j for j in range(n - c)],
                        ind="            ")
            o("            break;")
        o("        }")
    else:
        c = r_const // 64
        b = r_const % 64
        if c >= 1 and b == 0:
            vt = ["x%d" % (j + c) for j in range(n - c)] + ["UWORD(1)"]
            wt = ["y%d" % (j + c) for j in range(n - c)]
            masked_step(o, n, c, vt, wt, ind="        ")
        else:
            o("        const int b = %d;" % b)
            o("")
            o("        lt = MPN_RIGHT_SHIFT_LOW(UWORD(1), x%d, b);"
              % (n - 1))
            vw_shifted(o, n, c, "        ")
            masked_step(o, n, c, ["v%d" % j for j in range(n - c)],
                        ["w%d" % j for j in range(n - c)],
                        ind="        ")
    o("    }")
    o("")
    o("#undef AP")
    o("")
    # residual t = Y / X (one division; X = 1.x with the units limb)
    o("    /* t = Y 2^(FLINT_BITS %d) / X */" % n)
    for j in range(n):
        o("    S[%d] = x%d;" % (j, j))
    o("    S[%d] = 1;" % n)
    o("    nz = %s;" % " | ".join("(y%d != 0)" % j for j in range(n)))
    o("    if (nz)")
    o("    {")
    for j in range(n):
        o("        nd[%d] = 0;" % j)
    for j in range(n):
        o("        nd[%d] = y%d;" % (n + j, j))
    o("        mpn_tdiv_qr(t, nd, 0, nd, %d, S, %d);" % (2 * n, n + 1))
    o("    }")
    o("    else")
    o("    {")
    for j in range(n):
        o("        t[%d] = 0;" % j)
    o("    }")
    o("")
    if opt:
        o("    %s(q, t);" % series_name)
    else:
        o("    fixed_atan_rs(q, t, %d);" % n)
    o("")
    # res = acc + atan(t)
    aa = ["q[%d]" % j for j in range(n)]
    add_into(o, "    ", "a", aa, n)
    for j in range(n):
        o("    res[%d] = a%d;" % (j, j))
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
