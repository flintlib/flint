#!/usr/bin/env python3
"""Body emitter for dev/tune_mp_real.py (not run directly): the
dense-table per-size kernels, straight-line code on named scalars
against the tables of src/mp_real/dense_tables.c (see
dev/gen_mp_real_dense_tables.py for the reductions and their bounds).

exp, n fraction limbs, L levels (reduced radius r = 7 L - 1):

    t = v - L_1[j_1], j_1 the top seven bits of v
    for k = 2 .. L:  j_k = floor(t 2^7k),  t <- t - G_k[j_k]
    y = exp(t) by the hand series for t < 2^-r (y = 1 + F, F < 1)
    for k = L .. 1:  F <- F + (1 + F) c_k,  c_k = j_k 2^-7k (one limb),
                     c_1 = C_1[j_1] / 2^64 (the last one with a carry)
    exp(v) = 2^s (1 + F),  s = [j_1 >= S_1]

The single-limb products truncate (each step at most two ulps low, the
low word of the lowest partial product dropped).  With the top level
(r = 63, see RS_TIERS below), a tenth step takes the whole top limb of t
after nine levels as its index j_t, subtracts log(1 + j_t 2^-64) and
multiplies by 1 + j_t 2^-64 in the reconstruction.
"""

ADD = {2: "add_ssaaaa", 3: "add_sssaaaaaa", 4: "add_ssssaaaaaaaa",
       5: "add_sssssaaaaaaaaaa", 6: "add_ssssssaaaaaaaaaaaa",
       7: "add_sssssssaaaaaaaaaaaaaa",
       8: "add_ssssssssaaaaaaaaaaaaaaaa"}
SUB = {2: "sub_ddmmss", 3: "sub_dddmmmsss", 4: "sub_ddddmmmmssss",
       5: "sub_dddddmmmmmsssss", 6: "sub_ddddddmmmmmmssssss",
       7: "sub_dddddddmmmmmmmsssssss",
       8: "sub_ddddddddmmmmmmmmssssssss"}
for _k in (9, 10, 11):
    ADD[_k] = "add_" + "s" * _k + "a" * (2 * _k)
    SUB[_k] = "sub_" + "d" * _k + "m" * _k + "s" * _k

NMIN, NMAX = 4, 10
BITS = 7

# The series that finishes a dense kernel: by default the hand series of
# the reduced radius r = 7 L - 1 (L levels of seven bits); two values of
# r select instead the hardcoded rectangular splitting of exp_rs.c and
# trig_rs.c through _mp_real_exp_rs, _mp_real_atan_rs and
# _mp_real_sin_cos_rs, whose tiers sit at 2^-32 and 2^-64:
#
#   r = 32: five levels (the remainder below 2^-34), the 2^-32 tier;
#   r = 63: nine levels and a top level on the rest of the top limb (the
#           remainder below 2^-64 but for rare cases, see
#           gen_mp_real_dense_tables.py; exp only), the 2^-64 tier.
RS_TIERS = (32, 63)



def tabref(fam, k, j, n):
    """the address of the top n limbs of entry j (a C expression) of
    level k of fam: the first two levels from the W = 21 tables of
    dense_rt_tables.c (shared with the run-time kernels), the others
    from dense_tables.c"""
    if k <= 2:
        name = {("exp", 1): "exp_l1", ("exp", 2): "exp_g2",
                ("log", 1): "log_m1", ("log", 2): "log_m2",
                ("atan", 1): "atan_a1", ("atan", 2): "atan_a2",
                ("sc", 1): "sin_cos_t1"}[(fam, k)]
        return ("_mp_real_dense_rt_%s + %s * MP_REAL_DENSE_RT_W"
                " + (MP_REAL_DENSE_RT_W - %d)" % (name, j, n))
    tab, mac = {"exp": ("exp_g", "EXP_G"), "log": ("log_m", "LOG_M"),
                "atan": ("atan_a", "ATAN_A")}[fam]
    return ("_mp_real_dense_%s + MP_REAL_DENSE_%s%d + %s * MP_REAL_DENSE_W"
            " + (MP_REAL_DENSE_W - %d)" % (tab, mac, k, j, n))

def dense_levels(r):
    """(L, top, rs) for the kernel parameter r"""
    if r == 32:
        return 5, False, True
    if r == 63:
        return 9, True, True
    if (r + 1) % BITS != 0:
        raise ValueError("dense kernels take r = 7 L - 1, 32 or 63")
    return (r + 1) // BITS, False, False


def emit_rs_series(o, func, n, name):
    """the series wrapper of the hardcoded tiers (see RS_TIERS)"""
    o("/* the series: the hardcoded rectangular splitting of the 2^-32 and")
    o("   2^-64 tiers (%s) */" % ("exp_rs.c" if func == "exp" else "trig_rs.c"))
    o("FLINT_FORCE_INLINE void")
    if func == "sin_cos":
        o("%s(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)" % name)
        o("{")
        o("    _mp_real_sin_cos_rs(ysin, ycos, NULL, x, %d);" % n)
    else:
        o("%s(nn_ptr res, nn_srcptr x)" % name)
        o("{")
        o("    _mp_real_%s_rs(res, NULL, x, %d);" % (func, n))
    o("}")
    o("")


def chain(o, mac, dst, a, b, indent="    "):
    """dst = a op b on limb lists given most significant first"""
    o("%s%s(%s," % (indent, mac, ", ".join(dst)))
    o("%s    %s," % (indent, ", ".join(a)))
    o("%s    %s);" % (indent, ", ".join(b)))


def gen_exp(n, L, fn_name, series_name, static, top=False):
    assert NMIN <= n <= NMAX and 2 <= L and BITS * L <= 64
    assert not top or L == 9
    out = []
    o = out.append
    t = ["t%d" % i for i in range(n - 1, -1, -1)]          # msl first
    F = ["f%d" % i for i in range(n - 1, -1, -1)]
    if top:
        o("/* exp, n = %d: dense tables, %d levels of %d bits and a top level"
          % (n, L, BITS))
        o("   (the remainder below 2^-64 but for rare cases) */")
    else:
        o("/* exp, n = %d: dense tables, %d levels of %d bits (the remainder"
          % (n, L, BITS))
        o("   below 2^-%d) */" % (BITS * L - 1))
    o(("static void" if static else "void"))
    o("%s(nn_ptr res, nn_srcptr x)" % fn_name)
    o("{")
    o("    ulong %s;" % ", ".join("t%d" % i for i in range(n)))
    o("    ulong %s;" % ", ".join("f%d" % i for i in range(n)))
    o("    ulong %s;" % ", ".join("h%d" % i for i in range(n)))
    o("    ulong %s;" % ", ".join("l%d" % i for i in range(1, n)))
    o("    ulong %s, c, e, s;" % ", ".join(["j%d" % k for k in range(1, L + 1)]
                                         + (["jt"] if top else [])))
    o("    ulong ts[%d], ys[%d];" % (n, n + 1))
    o("    const ulong * T;")
    o("")
    o("    /* t = v - L_1[j_1], then t <- t - log(1 + j_k 2^-7k) */")
    o("    j1 = x[%d] >> (FLINT_BITS - %d);" % (n - 1, BITS))
    o("    T = %s;" % tabref("exp", 1, "j1", n))
    chain(o, SUB[n], t, ["x[%d]" % i for i in range(n - 1, -1, -1)],
          ["T[%d]" % i for i in range(n - 1, -1, -1)])
    for k in range(2, L + 1):
        o("    j%d = t%d >> (FLINT_BITS - %d);" % (k, n - 1, BITS * k))
        o("    T = %s;" % tabref("exp", k, "j%d" % k, n))
        chain(o, SUB[n], t, t, ["T[%d]" % i for i in range(n - 1, -1, -1)])
    if top:
        o("    /* the top level: t <- t - log(1 + j_t 2^-64) */")
        o("    jt = t%d;" % (n - 1))
        o("    T = _mp_real_dense_exp_gtop + jt * MP_REAL_DENSE_W"
          " + (MP_REAL_DENSE_W - %d);" % n)
        chain(o, SUB[n], t, t, ["T[%d]" % i for i in range(n - 1, -1, -1)])
    o("")
    o("    /* 1 + F = exp(t) */")
    for i in range(n):
        o("    ts[%d] = t%d;" % (i, i))
    o("    %s(ys, ts);" % series_name)
    for i in range(n):
        o("    f%d = ys[%d];" % (i, i))
    o("")
    o("    /* 1 + F <- (1 + F)(1 + c), c = j_k 2^-7k: F + c + F c, F c"
      " truncated */")
    for k in ([None] if top else []) + list(range(L, 1, -1)):
        if k is None:
            o("    c = jt;")
        else:
            o("    c = j%d << (FLINT_BITS - %d);" % (k, BITS * k))
        o("    h0 = n_mulhi(f0, c);")
        for i in range(1, n):
            o("    umul_ppmm(h%d, l%d, f%d, c);" % (i, i, i))
        chain(o, ADD[n], F, F, ["h%d" % i for i in range(n - 1, -1, -1)])
        chain(o, ADD[n], F, F,
              ["c"] + ["l%d" % i for i in range(n - 1, 0, -1)])
    o("")
    o("    /* (e, F) = (1 + F)(1 + C_1[j_1]) < 2.02, then times 2^s */")
    o("    c = _mp_real_dense_exp_c1[j1];")
    o("    h0 = n_mulhi(f0, c);")
    for i in range(1, n):
        o("    umul_ppmm(h%d, l%d, f%d, c);" % (i, i, i))
    chain(o, ADD[n + 1], ["e"] + F, ["UWORD(1)"] + F,
          ["UWORD(0)"] + ["h%d" % i for i in range(n - 1, -1, -1)])
    chain(o, ADD[n + 1], ["e"] + F, ["e"] + F,
          ["UWORD(0)", "c"] + ["l%d" % i for i in range(n - 1, 0, -1)])
    o("    s = (j1 >= MP_REAL_DENSE_EXP_S1);")
    o("    res[%d] = (e << s) | ((f%d >> 1) >> (FLINT_BITS - 1 - s));"
      % (n, n - 1))
    for i in range(n - 1, 0, -1):
        o("    res[%d] = (f%d << s) | ((f%d >> 1) >> (FLINT_BITS - 1 - s));"
          % (i, i, i - 1))
    o("    res[0] = f0 << s;")
    o("}")
    o("")
    return out


# ---------------------------------------------------------------------
# log1p
# ---------------------------------------------------------------------

def _trig_helpers():
    import importlib.util
    import os
    import sys
    name = "g_trig_series_helpers"
    if name in sys.modules:
        return sys.modules[name]
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "gen_mp_real_trig_hand_mixed.py")
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def log1p_series_levels(n, r):
    """V_0 .. V_{N-1}: log1p(x) = x - x^2 V_0, V_k = 1/(k+2) - x V_{k+1};
    the first dropped term x^(N+2)/(N+2) lies below 2^-64n"""
    from math import log2
    k = 0
    while -r * (k + 2) - log2(k + 2) >= -64 * n:
        k += 1
    return max(1, k)


def log1p_series_width(n, r, k):
    bits = 64 * n - 6 - r * (k + 2)
    if bits <= 0:
        return 1
    return max(1, -(-bits // 64))


def emit_log1p_series(o, n, r, name):
    """(res, n) = log(1 + x) for x < 2^-r at n limbs: Horner in x with
    the tapered widths of the hand series (an error in V_k reaches the
    result times x^(k+2)); every V_k is positive (x V_{k+1} < 1/(k+2))"""
    from fractions import Fraction
    h = _trig_helpers()
    N = log1p_series_levels(n, r)
    L = [min(n, log1p_series_width(n, r, k)) for k in range(N)]
    for k in range(N - 1, 0, -1):
        L[k - 1] = max(L[k - 1], L[k])
    o("/* log1p, n = %d, r = %d: %d levels, widths %s */" % (n, r, N, L))
    o("static void")
    o("%s(nn_ptr res, nn_srcptr x)" % name)
    o("{")
    o("    ulong z[%d], V[%d], T[%d], P[%d], Z[2];" % (n, n, n, n))
    o("    ulong v, w1, w0, h, l;")
    o("")
    o("    (void) v; (void) w1; (void) w0; (void) h; (void) l;")
    o("    (void) Z; (void) V; (void) T; (void) P;")
    o("")
    k = N - 1
    cur = L[k]
    c = h.cbits(Fraction(1, k + 2), cur)
    if cur == 1:
        o("    v = UWORD(0x%016x);" % c[0])
    elif cur == 2:
        o("    w1 = UWORD(0x%016x); w0 = UWORD(0x%016x);" % (c[1], c[0]))
    else:
        for i in range(cur):
            o("    V[%d] = UWORD(0x%016x);" % (i, c[i]))
    for k in range(N - 2, -1, -1):
        prev, cur = cur, L[k]
        c = h.cbits(Fraction(1, k + 2), cur)
        o("")
        o("    /* V_%d = 1/%d - x V_%d   (%d x %d limbs) */"
          % (k, k + 2, k + 1, cur, prev))
        if cur == 1:
            o("    v = UWORD(0x%016x) - n_mulhi(x[%d], v);" % (c[0], n - 1))
        elif cur == 2 and prev == 1:
            o("    _mp_real_mulhi_2x1(&h, &l, x[%d], x[%d], v);"
              % (n - 1, n - 2))
            o("    sub_ddmmss(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (c[1], c[0]))
        elif cur == 2:
            o("    _mp_real_mulhi_2x2_sloppy(&h, &l, x[%d], x[%d], w1, w0);"
              % (n - 1, n - 2))
            o("    sub_ddmmss(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (c[1], c[0]))
        else:
            if prev == 1:
                for i in range(cur - 1):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = v;" % (cur - 1))
            elif prev == 2:
                for i in range(cur - 2):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = w0; P[%d] = w1;" % (cur - 2, cur - 1))
            else:
                for i in range(cur - prev):
                    o("    P[%d] = UWORD(0);" % i)
                for i in range(prev):
                    o("    P[%d] = V[%d];" % (cur - prev + i, i))
            h.emit_mulhigh(o, "T", "x + %d" % (n - cur), "P", cur)
            for i in range(cur):
                o("    P[%d] = UWORD(0x%016x);" % (i, c[i]))
            h.emit_addsub(o, True, "V", "P", "T", cur)
    o("")
    o("    /* log1p = x - x^2 V_0 */")
    if cur == 1:
        for i in range(n - 1):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = v;" % (n - 1))
    elif cur == 2:
        for i in range(n - 2):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = w0; P[%d] = w1;" % (n - 2, n - 1))
    else:
        for i in range(n - cur):
            o("    P[%d] = UWORD(0);" % i)
        for i in range(cur):
            o("    P[%d] = V[%d];" % (n - cur + i, i))
    h.emit_sqrhigh(o, "z", "x", n)
    h.emit_mulhigh(o, "T", "z", "P", n)
    h.emit_addsub(o, True, "res", "x", "T", n)
    o("}")
    o("")


def log_d_expr(k):
    """C expression for the one-limb D_k(j) (k >= 2) from the index jk"""
    s = BITS * k
    if 2 * s <= 64:
        c2 = "((j%d * j%d) << %d)" % (k, k, 64 - 2 * s)
    else:
        sh = 2 * s - 64
        c2 = "((j%d * j%d + UWORD(%d)) >> %d)" % (k, k, (1 << sh) - 1, sh)
    return "(j%d << %d) - %s" % (k, 64 - s, c2)


def gen_log1p(n, L, fn_name, series_name, static):
    assert NMIN <= n <= NMAX and 2 <= L and BITS * L <= 64
    out = []
    o = out.append
    u = ["u%d" % i for i in range(n - 1, -1, -1)]
    a = ["a%d" % i for i in range(n - 1, -1, -1)]
    o("/* log1p, n = %d: dense tables, %d levels of %d bits (the remainder"
      % (n, L, BITS))
    o("   below 2^-%d) */" % (BITS * L - 1))
    o(("static void" if static else "void"))
    o("%s(nn_ptr res, nn_srcptr x)" % fn_name)
    o("{")
    o("    ulong %s;" % ", ".join("u%d" % i for i in range(n)))
    o("    ulong %s;" % ", ".join("a%d" % i for i in range(n)))
    o("    ulong %s;" % ", ".join("h%d" % i for i in range(n)))
    o("    ulong %s;" % ", ".join("l%d" % i for i in range(1, n)))
    o("    ulong %s, d;" % ", ".join("j%d" % k for k in range(1, L + 1)))
    o("    ulong ts[%d], ys[%d];" % (n, n))
    o("    const ulong * T;")
    o("")
    for i in range(n):
        o("    u%d = x[%d];" % (i, i))
    o("")
    for k in range(1, L + 1):
        o("    /* level %d: u <- u - D - u D (u D truncated), acc += M_%d */"
          % (k, k))
        o("    j%d = u%d >> (FLINT_BITS - %d);" % (k, n - 1, BITS * k))
        if k == 1:
            o("    d = _mp_real_dense_log_d1[j1];")
        else:
            o("    d = %s;" % log_d_expr(k))
        o("    h0 = n_mulhi(u0, d);")
        for i in range(1, n):
            o("    umul_ppmm(h%d, l%d, u%d, d);" % (i, i, i))
        chain(o, SUB[n], u, u, ["h%d" % i for i in range(n - 1, -1, -1)])
        chain(o, SUB[n], u, u, ["d"] + ["l%d" % i for i in range(n - 1, 0, -1)])
        o("    T = %s;" % tabref("log", k, "j%d" % k, n))
        if k == 1:
            for i in range(n):
                o("    a%d = T[%d];" % (i, i))
        else:
            chain(o, ADD[n], a, a, ["T[%d]" % i for i in range(n - 1, -1, -1)])
        o("")
    o("    /* + log(1 + u) */")
    for i in range(n):
        o("    ts[%d] = u%d;" % (i, i))
    o("    %s(ys, ts);" % series_name)
    chain(o, ADD[n], ["res[%d]" % i for i in range(n - 1, -1, -1)], a,
          ["ys[%d]" % i for i in range(n - 1, -1, -1)])
    o("}")
    o("")
    return out


# ---------------------------------------------------------------------
# atan
# ---------------------------------------------------------------------

def gen_atan(n, L, fn_name, series_name, static):
    """atan(v), v in [0, 1) at n limbs: rotations of (X, Y) = (1, v) by
    1 - i c_k (see dev/gen_mp_real_dense_tables.py), X = 1 + (x_i), then
    atan(Y/X) by one division and the series"""
    assert NMIN <= n <= NMAX and 2 <= L and BITS * L <= 64 - 9 + BITS
    out = []
    o = out.append
    X = ["x%d" % i for i in range(n - 1, -1, -1)]
    Y = ["y%d" % i for i in range(n - 1, -1, -1)]
    A = ["a%d" % i for i in range(n - 1, -1, -1)]
    o("/* atan, n = %d: dense tables, %d levels of %d bits (the remaining"
      % (n, L, BITS))
    o("   angle below 2^-%d) */" % (BITS * L - 1))
    o(("static void" if static else "void"))
    o("%s(nn_ptr res, nn_srcptr v)" % fn_name)
    o("{")
    for p in ("x", "y", "a", "hx", "hy"):
        o("    ulong %s;" % ", ".join("%s%d" % (p, i) for i in range(n)))
    for p in ("lx", "ly"):
        o("    ulong %s;" % ", ".join("%s%d" % (p, i) for i in range(1, n)))
    o("    ulong %s, c, R, rr, q;" % ", ".join("j%d" % k for k in range(1, L + 1)))
    o("    ulong S[%d], nd[%d], t[%d], ys[%d];" % (n + 1, 2 * n, n, n))
    o("    const ulong * T;")
    o("")
    o("    /* level 1: c = j_1/2^7 <= v, X = 1 + c v, Y = v - c */")
    o("    j1 = v[%d] >> (FLINT_BITS - %d);" % (n - 1, BITS))
    o("    c = j1 << (FLINT_BITS - %d);" % BITS)
    for i in range(n - 1):
        o("    y%d = v[%d];" % (i, i))
    o("    y%d = v[%d] & ((UWORD(1) << (FLINT_BITS - %d)) - 1);"
      % (n - 1, n - 1, BITS))
    o("    hx0 = n_mulhi(v[0], c);")
    for i in range(1, n):
        o("    umul_ppmm(hx%d, lx%d, v[%d], c);" % (i, i, i))
    chain(o, ADD[n], X, ["hx%d" % i for i in range(n - 1, -1, -1)],
          ["UWORD(0)"] + ["lx%d" % i for i in range(n - 1, 0, -1)])
    o("    T = %s;" % tabref("atan", 1, "j1", n))
    for i in range(n):
        o("    a%d = T[%d];" % (i, i))
    o("")
    o("    /* R ~ 2^64 / X (rounded so that R / 2^64 >= 1/X up to 2^-62) */")
    o("    udiv_qrnnd(R, rr, (UWORD(1) << (FLINT_BITS - 1)) - 1, ~UWORD(0),")
    o("        (UWORD(1) << (FLINT_BITS - 1)) | (x%d >> 1));" % (n - 1))
    o("    (void) rr;")
    for k in range(2, L + 1):
        s = BITS * k - 9
        o("")
        o("    /* level %d: j = floor(2^%d Y/X) from the estimate (lowered"
          " by 2^-12), then" % (k, BITS * k))
        o("       X <- X + c Y, Y <- Y - c X on the old values,"
          " c = j 2^-%d */" % (BITS * k))
        o("    q = n_mulhi((y%d << %d) | (y%d >> %d), R);"
          % (n - 1, s, n - 2, 64 - s))
        o("    q -= q >> 12;")
        o("    j%d = q >> %d;" % (k, 64 - 9))
        o("    c = j%d << (FLINT_BITS - %d);" % (k, BITS * k))
        o("    hx0 = n_mulhi(x0, c);")
        o("    hy0 = n_mulhi(y0, c);")
        for i in range(1, n):
            o("    umul_ppmm(hx%d, lx%d, x%d, c);" % (i, i, i))
            o("    umul_ppmm(hy%d, ly%d, y%d, c);" % (i, i, i))
        chain(o, ADD[n], X, X, ["hy%d" % i for i in range(n - 1, -1, -1)])
        chain(o, ADD[n], X, X,
              ["UWORD(0)"] + ["ly%d" % i for i in range(n - 1, 0, -1)])
        chain(o, SUB[n], Y, Y, ["hx%d" % i for i in range(n - 1, -1, -1)])
        chain(o, SUB[n], Y, Y, ["c"] + ["lx%d" % i for i in range(n - 1, 0, -1)])
        o("    T = %s;" % tabref("atan", k, "j%d" % k, n))
        chain(o, ADD[n], A, A, ["T[%d]" % i for i in range(n - 1, -1, -1)])
    o("")
    o("    /* t = Y/X < 2^-%d, + atan(t) */" % (BITS * L - 1))
    o("    /* 1/X to the precision t needs (t < 2^-%d), then t = Y (1/X) */"
      % (BITS * L - 1))
    for i in range(n):
        o("    S[%d] = x%d;" % (i, i))
    o("    _flint_mpn_recip_small(nd, S, %d, %d);" % (n, 64 * n - (BITS * L - 1) + 4))
    for i in range(n):
        o("    S[%d] = y%d;" % (i, i))
    o("    flint_mpn_mulhigh_n(t, S, nd, %d);" % n)
    o("    %s(ys, t);" % series_name)
    chain(o, ADD[n], ["res[%d]" % i for i in range(n - 1, -1, -1)], A,
          ["ys[%d]" % i for i in range(n - 1, -1, -1)])
    o("}")
    o("")
    return out


# ---------------------------------------------------------------------
# sin, cos
# ---------------------------------------------------------------------

def gen_sin_cos(n, L, fn_name, series_name, static):
    """(sin v, cos v), v in [0, 1) at n limbs (both n + 1 limbs, the
    top one the integral limb): the reduction of w = v/2 (see
    dev/gen_mp_real_dense_tables.py), (cos t, sin t) by the series, the
    rotations by 1 + i c_k, then the double-angle formulas divided by
    the squared norm N = X^2 + Y^2 (one reciprocal)"""
    assert NMIN <= n <= NMAX and 2 <= L and BITS * L <= 64
    out = []
    o = out.append
    t = ["t%d" % i for i in range(n - 1, -1, -1)]
    X = ["x%d" % i for i in range(n, -1, -1)]          # with integral limb
    Y = ["y%d" % i for i in range(n - 1, -1, -1)]
    def rev(p, a=n - 1, b=0):
        return ["%s%d" % (p, i) for i in range(a, b - 1, -1)]
    o("/* sin and cos, n = %d: dense tables, %d levels of %d bits on v/2 (the"
      % (n, L, BITS))
    o("   remainder below 2^-%d) */" % (BITS * L - 1))
    o(("static void" if static else "void"))
    o("%s(nn_ptr ysin, nn_ptr ycos, nn_srcptr v)" % fn_name)
    o("{")
    o("    ulong %s;" % ", ".join("t%d" % i for i in range(n)))
    o("    ulong %s;" % ", ".join("x%d" % i for i in range(n + 1)))
    for p in ("y", "hx", "hy"):
        o("    ulong %s;" % ", ".join("%s%d" % (p, i) for i in range(n)))
    for p in ("lx", "ly"):
        o("    ulong %s;" % ", ".join("%s%d" % (p, i) for i in range(1, n)))
    o("    ulong %s, c, mx, ni, di, si;" % ", ".join("j%d" % k for k in range(1, L + 1)))
    o("    ulong ts[%d], ss[%d], cc[%d], XF[%d], YF[%d], X2[%d], Y2[%d],"
      " XY[%d], NF[%d], I[%d];" % ((n, n + 1, n + 1) + (n,) * 7))
    o("    const ulong * T;")
    o("")
    o("    /* t = w = v/2 (the lowest bit of v dropped) */")
    for i in range(n - 1):
        o("    t%d = (v[%d] >> 1) | (v[%d] << (FLINT_BITS - 1));" % (i, i, i + 1))
    o("    t%d = v[%d] >> 1;" % (n - 1, n - 1))
    o("")
    o("    /* t <- t - atan(c_k), c_1 = C_1[j_1], c_k = j_k 2^-7k */")
    o("    j1 = t%d >> (FLINT_BITS - %d);" % (n - 1, BITS))
    o("    T = %s;" % tabref("sc", 1, "j1", n))
    chain(o, SUB[n], t, t, ["T[%d]" % i for i in range(n - 1, -1, -1)])
    for k in range(2, L + 1):
        o("    j%d = t%d >> (FLINT_BITS - %d);" % (k, n - 1, BITS * k))
        o("    T = %s;" % tabref("atan", k, "j%d" % k, n))
        chain(o, SUB[n], t, t, ["T[%d]" % i for i in range(n - 1, -1, -1)])
    o("")
    o("    /* (X, Y) = (cos t, sin t) */")
    for i in range(n):
        o("    ts[%d] = t%d;" % (i, i))
    o("    %s(ss, cc, ts);" % series_name)
    for i in range(n + 1):
        o("    x%d = cc[%d];" % (i, i))
    for i in range(n):
        o("    y%d = ss[%d];" % (i, i))
    o("")
    o("    /* (X, Y) <- (X - c Y, Y + c X), X = x_n + (fraction), x_n in {0, 1} */")
    for k in range(L, 0, -1):
        if k == 1:
            o("    c = _mp_real_dense_sin_cos_c1[j1];")
        else:
            o("    c = j%d << (FLINT_BITS - %d);" % (k, BITS * k))
        o("    mx = -x%d;" % n)
        o("    hx0 = n_mulhi(x0, c);")
        o("    hy0 = n_mulhi(y0, c);")
        for i in range(1, n):
            o("    umul_ppmm(hx%d, lx%d, x%d, c);" % (i, i, i))
            o("    umul_ppmm(hy%d, ly%d, y%d, c);" % (i, i, i))
        chain(o, SUB[n + 1], X, X, ["UWORD(0)"] + rev("hy"))
        chain(o, SUB[n + 1], X, X, ["UWORD(0)", "UWORD(0)"] + rev("ly", n - 1, 1))
        chain(o, ADD[n], Y, Y, rev("hx"))
        chain(o, ADD[n], Y, Y, ["(c & mx)"] + rev("lx", n - 1, 1))
    o("")
    o("    /* X^2, Y^2, XY, N = X^2 + Y^2 (in [1, 1.3] up to rounding) */")
    o("    mx = -x%d;" % n)
    for i in range(n):
        o("    XF[%d] = x%d;" % (i, i))
        o("    YF[%d] = y%d;" % (i, i))
    o("    flint_mpn_sqrhigh(X2, XF, %d);" % n)
    o("    flint_mpn_sqrhigh(Y2, YF, %d);" % n)
    o("    flint_mpn_mulhigh_n(XY, XF, YF, %d);" % n)
    # the masked operands of the chains below in memory (XF, YF are
    # free now): operands computed on the spot would need registers of
    # their own beside the chain's outputs, more than there are at low
    # optimization levels ("asm operand has impossible constraints")
    for i in range(n):
        o("    XF[%d] &= mx;" % i)
        o("    YF[%d] &= mx;" % i)
    # X^2 = x_n + 2 x_n F + F^2 -> (di, X2)
    o("    di = x%d;" % n)
    chain(o, ADD[n + 1], ["di"] + ["X2[%d]" % i for i in range(n - 1, -1, -1)],
          ["di"] + ["X2[%d]" % i for i in range(n - 1, -1, -1)],
          ["UWORD(0)"] + ["XF[%d]" % i for i in range(n - 1, -1, -1)])
    chain(o, ADD[n + 1], ["di"] + ["X2[%d]" % i for i in range(n - 1, -1, -1)],
          ["di"] + ["X2[%d]" % i for i in range(n - 1, -1, -1)],
          ["UWORD(0)"] + ["XF[%d]" % i for i in range(n - 1, -1, -1)])
    # XY += x_n Y (below 2/3, no integral part); 2XY -> (si, XY)
    chain(o, ADD[n], ["XY[%d]" % i for i in range(n - 1, -1, -1)],
          ["XY[%d]" % i for i in range(n - 1, -1, -1)],
          ["YF[%d]" % i for i in range(n - 1, -1, -1)])
    o("    si = XY[%d] >> (FLINT_BITS - 1);" % (n - 1))
    o("    mpn_lshift(XY, XY, %d, 1);" % n)
    # N = X^2 + Y^2 -> (ni, NF); D = X^2 - Y^2 -> (di, X2)
    chain(o, ADD[n + 1], ["ni"] + ["NF[%d]" % i for i in range(n - 1, -1, -1)],
          ["di"] + ["X2[%d]" % i for i in range(n - 1, -1, -1)],
          ["UWORD(0)"] + ["Y2[%d]" % i for i in range(n - 1, -1, -1)])
    chain(o, SUB[n + 1], ["di"] + ["X2[%d]" % i for i in range(n - 1, -1, -1)],
          ["di"] + ["X2[%d]" % i for i in range(n - 1, -1, -1)],
          ["UWORD(0)"] + ["Y2[%d]" % i for i in range(n - 1, -1, -1)])
    o("    /* N >= 1 exactly; rounded below 1, take N = 1 */")
    o("    mx = -ni;")
    for i in range(n):
        o("    NF[%d] &= mx;" % i)
    o("    _flint_mpn_recip_small(I, NF, %d, %d);" % (n, 64 * n - 3))
    o("")
    o("    /* cos v = D/N, sin v = 2XY/N (D, 2XY with integral limbs di, si) */")
    o("    flint_mpn_mulhigh_n(ycos, X2, I, %d);" % n)
    o("    flint_mpn_mulhigh_n(ysin, XY, I, %d);" % n)
    o("    mx = -di;")
    for i in range(n):
        o("    XF[%d] = I[%d] & mx;" % (i, i))
    o("    mx = -si;")
    for i in range(n):
        o("    YF[%d] = I[%d] & mx;" % (i, i))
    chain(o, ADD[n + 1], ["ycos[%d]" % i for i in range(n, -1, -1)],
          ["UWORD(0)"] + ["ycos[%d]" % i for i in range(n - 1, -1, -1)],
          ["UWORD(0)"] + ["XF[%d]" % i for i in range(n - 1, -1, -1)])
    chain(o, ADD[n + 1], ["ysin[%d]" % i for i in range(n, -1, -1)],
          ["UWORD(0)"] + ["ysin[%d]" % i for i in range(n - 1, -1, -1)],
          ["UWORD(0)"] + ["YF[%d]" % i for i in range(n - 1, -1, -1)])
    o("}")
    o("")
    return out
