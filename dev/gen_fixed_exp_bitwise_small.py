#!/usr/bin/env python3
"""Regenerate src/fixed/exp_bitwise_rs_small.inc (run from the top-level
FLINT directory:
python3 dev/gen_fixed_exp_bitwise_small.py > src/fixed/exp_bitwise_rs_small.inc):
per-size register implementations
of the bitwise-reduction exp for n = 1..5, using the fixed-size
add/sub chains from longlong.h and MPN_RIGHT_SHIFT_LOW from
mpn_extras.h.  See exp_bitwise_rs.c for the algorithm."""

SUBN = {2: "sub_ddmmss", 3: "sub_dddmmmsss", 4: "sub_ddddmmmmssss",
        5: "sub_dddddmmmmmsssss", 6: "sub_ddddddmmmmmmssssss",
        7: "sub_dddddddmmmmmmmsssssss"}
SUB = {2: "sub_ddmmss", 3: "sub_dddmmmsss", 4: "sub_ddddmmmmssss",
       5: "sub_dddddmmmmmsssss", 6: "sub_ddddddmmmmmmssssss",
       7: "sub_dddddddmmmmmmmsssssss",
       8: "sub_ddddddddmmmmmmmmssssssss"}
ADD = {2: "add_ssaaaa", 3: "add_sssaaaaaa", 4: "add_ssssaaaaaaaa",
       5: "add_sssssaaaaaaaaaa", 6: "add_ssssssaaaaaaaaaaaa",
       7: "add_sssssssaaaaaaaaaaaaaa",
       8: "add_ssssssssaaaaaaaaaaaaaaaa"}


def lst(p, k, hi_to_lo=True):
    r = range(k - 1, -1, -1) if hi_to_lo else range(k)
    return ", ".join("%s%d" % (p, i) for i in r)


# Fully specialized sizes: the reduction parameter is a compile-time
# constant (so the clamps fold away and every window loop gets a fixed
# trip count) and the series is the matching one from
# dev/gen_fixed_exp_rs_hard.py (WINNERS_OPT), whose term count spends
# the 1/N! of the last coefficient instead of banking it.
#
# Note that a constant r only pays while r is SMALL.  Specializing at
# r = 32 measured markedly slower than passing 32 at run time (194 vs
# 141 ns at n = 4): the compiler fully unrolls the 32-iteration window,
# and the resulting code costs more than the folding saves.  These
# entries all sit at r <= 16, where unrolling is still a win.
OPT_R = {1: 12, 2: 16, 3: 16, 4: 16, 5: 16}


def gen(n, opt=False, fn_name=None, series_name=None, r_const=None,
        static=True):
    """opt=True with fn_name/series_name/r_const emits a fully
    specialized body under an arbitrary name with an arbitrary series,
    for the per-size production files and the tuner's candidates."""
    wn = n
    L = []
    o = L.append
    if fn_name is None:
        fn_name = ("_fixed_exp_bitwise_rs_opt_%d" % n) if opt \
            else ("_fixed_exp_bitwise_rs_%d" % n)
    if series_name is None and opt:
        series_name = "_fixed_exp_rs_opt_%d" % n
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
    o("    slong i, j, num, bj, nc;")
    o("    ulong %s, bw, h, e;" % lst("t", wn, False))
    o("    ulong %s;" % lst("d", wn, False))
    o("    ulong %s;" % lst("y", wn + 1, False))
    o("    ulong %s;" % lst("s", wn + 1, False))
    o("    ulong ts[%d], ys[%d], sh[%d];" % (wn, wn + 2, wn + 1))
    o("    slong used[%d * FLINT_BITS + 2];" % n)
    o("")
    # r < 32 shortens the reduction at the price of a longer series;
    # the wider-range series (x < 2^-16) exists for n <= 5.  In the
    # specialized variants r is already a compile-time constant, so
    # the clamps are omitted entirely.
    if not opt:
        if n <= 5:
            o("    r = FLINT_MAX(r, 16);")
        else:
            o("    r = FLINT_MAX(r, 32);")
        o("    r = FLINT_MIN((slong) r, FLINT_BITS * %d - 16);" % n)
    o("")
    o("    _fixed_exp_logs_ensure(%d, r);" % wn)
    o("    nc = _fixed_exp_logs_n;")
    o("")
    for k in range(wn):
        o("    t%d = x[%d];" % (k, k))
    o("")
    o("/* exact compare-subtract of L_ii, recording the index */")
    o("#define STEP(ii) \\")
    o("    do { \\")
    o("        nn_srcptr Lq = _fixed_exp_logs + (ii) * nc + (nc - %d); \\"
      % wn)
    o("        %s(bw, %s, \\" % (SUB[wn + 1], lst("d", wn)))
    o("            UWORD(0), %s, \\" % lst("t", wn))
    o("            UWORD(0), %s); \\"
      % ", ".join("Lq[%d]" % k for k in range(wn - 1, -1, -1)))
    for k in range(wn):
        o("        t%d = (t%d & bw) | (d%d & ~bw); \\" % (k, k, k))
    o("        used[num] = (ii); \\")
    o("        num += (bw == 0); \\")
    o("    } while (0)")
    o("")
    o("/* apply the pending recorded subtractions */")
    o("#define FLUSH() \\")
    o("    for (; bj < num; bj++) \\")
    o("    { \\")
    o("        nn_srcptr Lq = _fixed_exp_logs + used[bj] * nc + (nc - %d); \\"
      % wn)
    if wn == 1:
        o("        t0 -= Lq[0]; \\")
    else:
        o("        %s(%s, \\" % (SUBN[wn], lst("t", wn)))
        o("            %s, \\" % lst("t", wn))
        o("            %s); \\"
      % ", ".join("Lq[%d]" % k for k in range(wn - 1, -1, -1)))
    o("    }")
    o("")
    o("    num = 0;")
    o("    bj = 0;")
    o("")
    if n <= 2:
        # at these widths the exact step is a two/three-limb chain and
        # the windowing overhead does not amortize: reduce directly
        o("    for (i = 0; i <= r; i++)")
        o("        STEP(i);")
        o("    bj = num;")
        o("    (void) h; (void) e;")
        o("")
    else:
        o("/* one reduction window: decisions on the single-limb model")
        o("   hreg, exact steps at the boundary and inside the ambiguity")
        o("   band (see _fixed_exp_reduce) */")
        o("#define WINDOW(cc, hreg) \\")
        o("    do { \\")
        if r_const is None:
            o("        if (FLINT_BITS * (cc) <= r) \\")
        o("        { \\")
        o("            slong i1 = FLINT_MIN((slong) r, \\")
        o("                FLINT_BITS * ((cc) + 1) - 1); \\")
        o("            nn_srcptr lp; \\")
        o("            FLUSH(); \\")
        o("            STEP(FLINT_BITS * (cc)); \\")
        o("            bj = num; \\")
        o("            h = hreg; \\")
        o("            e = 0; \\")
        o("            lp = _fixed_exp_logs + (FLINT_BITS * (cc) + 1) * nc \\")
        o("                + (nc - 1 - (cc)); \\")
        o("            for (i = FLINT_BITS * (cc) + 1; i <= i1; i++, lp += nc) \\")
        o("            { \\")
        o("                ulong lt = lp[0]; \\")
        o("                if (h - lt + e <= 2 * e) \\")
        o("                { \\")
        o("                    FLUSH(); \\")
        o("                    STEP(i); \\")
        o("                    bj = num; \\")
        o("                    h = hreg; \\")
        o("                    e = 0; \\")
        o("                } \\")
        o("                else \\")
        o("                { \\")
        o("                    ulong sub = (h > lt); \\")
        o("                    h -= (lt + 1) & (-sub); \\")
        o("                    used[num] = i; \\")
        o("                    num += sub; \\")
        o("                    e += sub; \\")
        o("                } \\")
        o("            } \\")
        o("        } \\")
        o("    } while (0)")
        o("")
        cmax = n - 1 if r_const is None else min(n - 1, r_const // 64)
        for c in range(0, cmax + 1):
            o("    WINDOW(%d, t%d);" % (c, wn - 1 - c))
        o("")
        o("#undef WINDOW")
        o("")
    o("    FLUSH();")
    o("")
    o("    /* one extra subtraction restores t < L_r < 2^-r */")
    o("    STEP(r);")
    o("")
    o("#undef STEP")
    o("#undef FLUSH")
    o("")
    o("    /* exp of the reduced argument */")
    for k in range(wn):
        o("    ts[%d] = t%d;" % (k, k))
    if opt:
        o("    %s(ys, ts);" % series_name)
    else:
        o("    fixed_exp_rs(ys, ts, %d);" % wn)
    for k in range(wn + 1):
        o("    y%d = ys[%d];" % (k, k))
    o("")
    o("    j = 0;")
    o("")
    o("    if (num > 0 && used[0] == 0)")
    o("    {")
    o("        %s(%s," % (ADD[wn + 1], lst("y", wn + 1)))
    o("            %s, %s);" % (lst("y", wn + 1), lst("y", wn + 1)))
    o("        j = 1;")
    o("    }")
    o("")
    o("    for (; j < num && used[j] < FLINT_BITS; j++)")
    o("    {")
    o("        int b = (int) used[j];")
    o("")
    for k in range(wn):
        o("        s%d = MPN_RIGHT_SHIFT_LOW(y%d, y%d, b);" % (k, k + 1, k))
    o("        s%d = y%d >> b;" % (wn, wn))
    o("        %s(%s," % (ADD[wn + 1], lst("y", wn + 1)))
    o("            %s, %s);" % (lst("y", wn + 1), lst("s", wn + 1)))
    o("    }")
    o("")
    for k in range(wn + 1):
        o("    ys[%d] = y%d;" % (k, k))
    o("")
    o("    /* rare factors with i >= FLINT_BITS */")
    o("    if (j < num)")
    o("        _fixed_exp_recon(ys, sh, %d, used, j, num);" % (wn + 1))
    o("")
    for k in range(wn + 1):
        o("    res[%d] = ys[%d];" % (k, k))
    o("}")
    o("")
    return "\n".join(L)


out = ["/* Generated by dev/gen_fixed_exp_bitwise_small.py -- do not edit.",
       "",
       "   Register implementations of fixed_exp_bitwise_rs for",
       "   n = 1..5; see exp_bitwise_rs.c for the algorithm. */", ""]
for n in range(1, 8):
    out.append(gen(n))
for n in sorted(OPT_R):
    out.append(gen(n, opt=True))
out.append("""/* sizes with a hardcoded reduction parameter (see OPT_R) */
static int
_fixed_exp_bitwise_rs_opt(nn_ptr res, nn_srcptr x, slong n)
{
    switch (n)
    {
        case 1: _fixed_exp_bitwise_rs_opt_1(res, x); return 1;
        case 2: _fixed_exp_bitwise_rs_opt_2(res, x); return 1;
        case 3: _fixed_exp_bitwise_rs_opt_3(res, x); return 1;
        case 4: _fixed_exp_bitwise_rs_opt_4(res, x); return 1;
        case 5: _fixed_exp_bitwise_rs_opt_5(res, x); return 1;
        default: return 0;
    }
}

static void
_fixed_exp_bitwise_rs_small(nn_ptr res, nn_srcptr x, slong n, int r)
{
    switch (n)
    {
        case 1: _fixed_exp_bitwise_rs_1(res, x, r); break;
        case 2: _fixed_exp_bitwise_rs_2(res, x, r); break;
        case 3: _fixed_exp_bitwise_rs_3(res, x, r); break;
        case 4: _fixed_exp_bitwise_rs_4(res, x, r); break;
        case 5: _fixed_exp_bitwise_rs_5(res, x, r); break;
        case 6: _fixed_exp_bitwise_rs_6(res, x, r); break;
        default: _fixed_exp_bitwise_rs_7(res, x, r); break;
    }
}
""")
print("\n".join(out).rstrip("\n"))


if __name__ == "__main__":
    import sys
    sys.exit("this module is a body-emitter library for "
             "dev/tune_fixed.py; run that instead")
