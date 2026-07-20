#!/usr/bin/env python3
"""Regenerate src/fixed/tan_halfangle_small.inc (run from the top-level
FLINT directory: python3 dev/gen_fixed_tan_halfangle_small.py >
src/fixed/tan_halfangle_small.inc).

The half-angle reconstruction tail for n = 2..4 with everything except
the reciprocal division in registers: the generic path spends a dozen
mpn calls (two rectangular products, adds, subtractions, shifts, the
output mulhighs) whose dispatch and buffer traffic dominate at these
sizes -- the same lesson the n = 1 scalar path and the register
rotation taught.  Here the products go through the hand 2x2 forms at
n = 2 and the assembly-backed flint_mpn_mulhigh_n/sqrhigh dispatch at
n = 3, 4; sums, differences and shifts are NN_ADD_k / NN_SUB_k chains
and unrolled funnel shifts; the one reciprocal per output stays an
mpn_tdiv_qr against a normalized divisor.

The algebra is the generic tail's, with the tangent t = TT/DE never
materialized:

    A = TT DE, B = DE^2, C = TT^2, S = B + C
    sin = 2A/S, cos = 1 - 2C/S, tan = 2A/(B - C)

Only r = _fixed_tan_opt_r[n] is supported (the r = 0 dispatch): the
series calls are the fixed-r tabulated ones.

The sloppy products cost a few ulp more than the generic path's exact
rectangular multiplies, well within the 6r + 128 / 8r + 256 budgets;
validate against MPFR after regenerating (the roundings differ from
the generic path BY DESIGN, so compare in ulp, not bits).
"""

NMIN, NMAX = 2, 4


def o_shift_right1(o, var, n, indent):
    """halve an n-limb register array in place"""
    pre = " " * indent
    for k in range(n - 1):
        o(pre + "%s[%d] = (%s[%d] >> 1) | (%s[%d] << (FLINT_BITS - 1));"
          % (var, k, var, k, var, k + 1))
    o(pre + "%s[%d] >>= 1;" % (var, n - 1))


def o_shift_left_var(o, dst, src, n, cnt, indent, unit=False):
    """dst[0..n-1(+unit)] = src[0..n-1] << cnt, cnt in [1, FLINT_BITS)"""
    pre = " " * indent
    if unit:
        o(pre + "%s[%d] = %s[%d] >> (FLINT_BITS - %s);" % (dst, n, src, n - 1, cnt))
    for k in range(n - 1, 0, -1):
        o(pre + "%s[%d] = (%s[%d] << %s) | (%s[%d] >> (FLINT_BITS - %s));"
          % (dst, k, src, k, cnt, src, k - 1, cnt))
    o(pre + "%s[0] = %s[0] << %s;" % (dst, src, cnt))


def emit_mulhi(o, res, a, b, n, indent):
    pre = " " * indent
    if n == 2:
        o(pre + "_fixed_mulhi_2x2_sloppy(%s + 1, %s, %s[1], %s[0], %s[1], %s[0]);"
          % (res, res, a, a, b, b))
    else:
        o(pre + "flint_mpn_mulhigh_n(%s, %s, %s, %d);" % (res, a, b, n))


def emit_sqrhi(o, res, a, n, indent):
    pre = " " * indent
    if n == 2:
        o(pre + "_fixed_sqrhi_2x2(%s + 1, %s, %s[1], %s[0]);" % (res, res, a, a))
    else:
        o(pre + "flint_mpn_sqrhigh(%s, %s, %d);" % (res, a, n))


def emit(o, n, r_const=None, fn_name=None, series_name=None,
         static=True):
    wn = n + 1
    if fn_name is None:
        fn_name = "_fixed_tan_halfangle_%d" % n
    if series_name is None:
        series_name = "_fixed_tan_rs_opt_%d" % n
    o("static void" if static else "void")
    o("%s(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan," % fn_name)
    o("    nn_srcptr x)")
    o("{")
    o("    slong num, nc;")
    if r_const is None:
        o("    int r = _fixed_tan_opt_r[%d], c;" % n)
    else:
        o("    const int r = %d;" % r_const)
        o("    int c;")
    o("    ulong t[%d], u[%d], wx[%d], wy[%d];" % (n, n, wn, wn))
    o("    ulong T[%d], D[%d], h[%d], A[%d], B[%d], C[%d];" % (n, wn, n, n, n, n))
    o("    ulong S[%d], N[%d], R[%d], rem[%d], v[%d];" % (n, 2 * n, wn, n, n))
    o("    slong * used;")
    o("    TMP_INIT;")
    o("")
    o("    _fixed_atans_ensure(%d, r);" % n)
    o("    nc = _fixed_atans_n;")
    o("")
    o("    TMP_START;")
    o("    used = TMP_ALLOC(FIXED_BITWISE_REDUCE_USED_ALLOC(r)"
      " * sizeof(slong));")
    o("")
    for k in range(n - 1):
        o("    t[%d] = (x[%d] >> 1) | (x[%d] << (FLINT_BITS - 1));" % (k, k, k + 1))
    o("    t[%d] = x[%d] >> 1;" % (n - 1, n - 1))
    o("    num = _fixed_bitwise_reduce(t, %d, r, 1, _fixed_atans, nc, used);" % n)
    o("    _fixed_tan_rotate_%d(wx, wy, used, num);" % n)
    o("    %s(u, t);          /* u = tan(t') */" % series_name)
    o("")
    o("    /* TT = wy + wx u < 0.64 (n limbs) */")
    emit_mulhi(o, "h", "wx", "u", n, 4)
    o("    if (wx[%d])" % n)
    o("        NN_ADD_%d(h, h, u);" % n)
    o("    NN_ADD_%d(T, h, wy);" % n)
    o("    /* DE = wx - wy u in (0.72, 1.17) */")
    emit_mulhi(o, "h", "wy", "u", n, 4)
    args_d = ", ".join("D[%d]" % k for k in range(n, -1, -1))
    args_x = ", ".join("wx[%d]" % k for k in range(n, -1, -1))
    args_h = "UWORD(0), " + ", ".join("h[%d]" % k for k in range(n - 1, -1, -1))
    o("    %s(%s," % (SUB[wn], args_d))
    o("        %s," % args_x)
    o("        %s);" % args_h)
    o("")
    o("    /* halved together when DE reaches one: DE normalized */")
    o("    if (D[%d])" % n)
    o("    {")
    o_shift_right1(o, "D", n, 8)
    o("        D[%d] |= UWORD(1) << (FLINT_BITS - 1);" % (n - 1))
    o_shift_right1(o, "T", n, 8)
    o("    }")
    o("")
    o("    /* A = TT DE, B = DE^2, C = TT^2 (A unconditionally: proving")
    o("       the sin/tan branches imply its initialization is beyond the")
    o("       compiler, and a cos-only caller wastes one cheap mulhigh) */")
    emit_mulhi(o, "A", "T", "D", n, 4)
    emit_sqrhi(o, "B", "D", n, 4)
    emit_sqrhi(o, "C", "T", n, 4)
    o("")
    for k in range(2 * n - 1):
        o("    N[%d] = ~UWORD(0);" % k)
    o("    N[%d] = ~UWORD(0) >> 1;" % (2 * n - 1))
    o("")
    o("    if (ysin != NULL || ycos != NULL)")
    o("    {")
    o("        /* S = B + C = DE^2 (1 + t^2), normalized to [0.5, 1):")
    o("           c = 2 - e is the doubling shift of sin = 2^(2-e) A R */")
    args_s = ", ".join("S[%d]" % k for k in range(n - 1, -1, -1))
    args_b = ", ".join("B[%d]" % k for k in range(n - 1, -1, -1))
    args_c = ", ".join("C[%d]" % k for k in range(n - 1, -1, -1))
    o("        ulong cy;")
    o("        %s(cy, %s," % (ADDC[n], args_s))
    o("            UWORD(0), %s," % args_b)
    o("            UWORD(0), %s);" % args_c)
    o("        if (cy)")
    o("        {")
    o_shift_right1(o, "S", n, 12)
    o("            S[%d] |= UWORD(1) << (FLINT_BITS - 1);" % (n - 1))
    o("            c = 1;")
    o("        }")
    o("        else if (!(S[%d] >> (FLINT_BITS - 1)))" % (n - 1))
    o("        {")
    o_shift_left_var(o, "S", "S", n, "1", 12)
    o("            c = 3;")
    o("        }")
    o("        else")
    o("            c = 2;")
    o("")
    o("        mpn_tdiv_qr(R, rem, 0, N, %d, S, %d);" % (2 * n, n))
    o("")
    o("        if (ysin != NULL)")
    o("        {")
    emit_mulhi(o, "v", "A", "R", n, 12)
    o_shift_left_var(o, "ysin", "v", n, "c", 12, unit=True)
    o("        }")
    o("        if (ycos != NULL)")
    o("        {")
    emit_mulhi(o, "v", "C", "R", n, 12)
    o_shift_left_var(o, "v", "v", n, "c", 12)
    o("            if (%s)" % " && ".join("v[%d] == 0" % k for k in range(n)))
    o("            {")
    for k in range(n):
        o("                ycos[%d] = 0;" % k)
    o("                ycos[%d] = 1;   /* cos exactly one */" % n)
    o("            }")
    o("            else")
    o("            {")
    args_yc = ", ".join("ycos[%d]" % k for k in range(n - 1, -1, -1))
    args_z = ", ".join("UWORD(0)" for k in range(n))
    args_v = ", ".join("v[%d]" % k for k in range(n - 1, -1, -1))
    o("                %s(%s," % (SUB[n], args_yc))
    o("                    %s," % args_z)
    o("                    %s);" % args_v)
    o("                ycos[%d] = 0;" % n)
    o("            }")
    o("        }")
    o("    }")
    o("")
    o("    if (ytan != NULL)")
    o("    {")
    o("        /* W = B - C = DE^2 (1 - t^2) > 0.175: at most two leading")
    o("           zero bits; tan = 2^(2+e) A R < 1.56 */")
    args_w = ", ".join("S[%d]" % k for k in range(n - 1, -1, -1))
    o("        %s(%s," % (SUB[n], args_w))
    o("            %s," % args_b)
    o("            %s);" % args_c)
    o("        c = flint_clz(S[%d]);" % (n - 1))
    o("        if (c)")
    o("        {")
    o_shift_left_var(o, "S", "S", n, "c", 12)
    o("        }")
    o("        c += 2;")
    o("")
    o("        mpn_tdiv_qr(R, rem, 0, N, %d, S, %d);" % (2 * n, n))
    o("")
    emit_mulhi(o, "v", "A", "R", n, 8)
    o_shift_left_var(o, "ytan", "v", n, "c", 8, unit=True)
    o("    }")
    o("")
    o("    TMP_END;")
    o("}")
    o("")


SUB = {2: "sub_ddmmss", 3: "sub_dddmmmsss", 4: "sub_ddddmmmmssss",
       5: "sub_dddddmmmmmsssss"}
ADDC = {2: "add_sssaaaaaa", 3: "add_ssssaaaaaaaa", 4: "add_sssssaaaaaaaaaa"}



def emit_named(n, r, fn_name, series_name, static):
    L = []
    emit(L.append, n, r_const=r, fn_name=fn_name,
         series_name=series_name, static=static)
    return L

if __name__ == "__main__":
    import sys
    sys.exit("this module is a body-emitter library for "
             "dev/tune_fixed.py; run that instead")
