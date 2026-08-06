#!/usr/bin/env python3
"""Regenerate src/fixed/trig_rs_hard.inc (run from the top-level
FLINT directory: python3 dev/gen_fixed_trig_rs_hard.py > src/fixed/trig_rs_hard.inc)."""

#!/usr/bin/env python3
"""
Generator for hardcoded fixed-point evaluators of the odd/even Taylor
families after argument reduction, for 0 <= x < 2^-zbits, zbits in
{64, 32}:

    void mpn_atan_series[32]_{xn}(mp_ptr res, mp_srcptr x)
    void mpn_atanh_series[32]_{xn}(mp_ptr res, mp_srcptr x)
    void mpn_sin_cos_series[32]_{xn}(mp_ptr ysin, mp_ptr ycos, mp_srcptr x)
    void mpn_sinh_cosh_series[32]_{xn}(mp_ptr ysinh, mp_ptr ycosh, mp_srcptr x)

atan/atanh: res[0..xn-1] ~= sum_k (-+1)^k x^(2k+1)/(2k+1).
sin/cos and sinh/cosh: ysin[0..xn] (top limb 0) and ycos[0..xn] get the
odd/even factorial sums; either output pointer may be NULL, in which
case all work specific to that series is skipped (the z powers are
shared).  x has xn fractional limbs; for zbits = 64 the top limb is
read by the final multiplications and must be zero.

Structure: substitute z = x^2 (power z^k carries Lz = 2*zbits/64
leading zero limbs exactly), rectangular splitting with even m and
descending plain mpn_addmul_1 / mpn_submul_1.  For the alternating
functions the accumulator may wrap negative between the dominant
units terms; the fresh gap and carry limbs then use the running
sign-extension limb fl_s, updated branch-free from each op's returned
borrow/carry (RS_* macros): a subtracted term can be exactly zero
(all z reads zero for very small x), so the sign cannot be tracked
statically.  All-positive stretches and the non-alternating functions
emit pure exp-style chains.  Even m keeps every boundary
multiplication and final division behind a dominant nonnegative value.

Term counts are per series: odd series take k while 2k+1 < xn
(zbits 64; k < xn for 32), even series while 2k < xn.  The k <= 1
even-series head is exact: ycos = 1 - z/2 by shift+negate,
ycosh = 1 + z/2, with the splitting only over k >= 2.  The
non-alternating series exclude k = 0 from the splitting entirely
(their sums exceed 1), finishing with res = x + mulhigh(q, x) and an
unconditional add.  The alternating odd series keep k = 0 as the
units value of the shared denominator.

Divisions by the shared denominator (3*5*...*33 for atan/atanh, 20!
for the factorial families) use either mpn_divrem_1 or a mulhigh by a
truncated inverse: the top n limbs of floor(2^(64*NMAX)/D) equal
floor(2^(64*n)/D) exactly (floor composition), so one constant array
serves every length, with a one-sided error <= 3 units.  The
truncated inverse strictly underestimates, so the quotient never
reaches 2^(64*xn); the divrem variant subtracts 1 from the
accumulator first (error < 2^-61 quotient units) for the same
guarantee, keeping all output assembly unconditional and branch-free.

Usage:
  python3 gen_atan_sincos.py all[32] > file.c
  python3 gen_atan_sincos.py best[32] WINNERS.json > file.c
"""

import sys
import json

FAC20 = 1
for i in range(1, 21):
    FAC20 *= i
D_ATAN = 1
for o in range(3, 35, 2):
    D_ATAN *= o
assert D_ATAN < 2 ** 64 and FAC20 < 2 ** 64

ATAN_NMAX = 17
SC_NMAX = 10
INVA_LEN = 36
INVF_LEN = 23


def limbs(v, n):
    return [(v >> (64 * i)) & (2 ** 64 - 1) for i in range(n)]


def emit_tables():
    out = []
    out.append("/* atan/atanh: c_k = D_A/(2k+1), D_A = 3*5*...*33 */")
    out.append("static const mp_limb_t atan_series_cs[%d] = {" % ATAN_NMAX)
    for k in range(ATAN_NMAX):
        out.append("    UWORD(%d)," % (D_ATAN // (2 * k + 1)))
    out.append("};")
    out.append("#define ATAN_SERIES_D UWORD(%d)" % D_ATAN)
    out.append("")
    out.append("/* sin/sinh: 20!/(2k+1)!, cos/cosh: 20!/(2k)! */")
    out.append("static const mp_limb_t sin_series_cs[%d] = {" % SC_NMAX)
    for k in range(SC_NMAX):
        c = FAC20
        for i in range(1, 2 * k + 2):
            c //= i
        out.append("    UWORD(%d)," % c)
    out.append("};")
    out.append("static const mp_limb_t cos_series_cs[%d] = {" % (SC_NMAX + 1))
    for k in range(SC_NMAX + 1):
        c = FAC20
        for i in range(1, 2 * k + 1):
            c //= i
        out.append("    UWORD(%d)," % c)
    out.append("};")
    out.append("#define SIN_COS_SERIES_D UWORD(%d)" % FAC20)
    out.append("")
    out.append("/* truncated inverses: the top n limbs equal")
    out.append("   floor(2^(64n)/D) exactly, for any n */")
    out.append("static const mp_limb_t atan_series_dinv[%d] = {" % INVA_LEN)
    for l in limbs((2 ** (64 * INVA_LEN)) // D_ATAN, INVA_LEN):
        out.append("    UWORD(%d)," % l)
    out.append("};")
    out.append("static const mp_limb_t sin_cos_series_dinv[%d] = {" % INVF_LEN)
    for l in limbs((2 ** (64 * INVF_LEN)) // FAC20, INVF_LEN):
        out.append("    UWORD(%d)," % l)
    out.append("};")
    out.append("""
/* one descending rectangular-splitting term over a possibly wrapped
   (two's complement) accumulator: op over a[0..n-1], fresh top limb
   a[n] absorbs the borrow/carry against the running sign-extension
   limb fl, updated branch-free */
#define RS_SUBMUL(a, n, src, c, fl) \\
    do { mp_limb_t __cy = mpn_submul_1(a, src, n, c); \\
         (a)[n] = (fl) - __cy; (fl) -= ((fl) < __cy); } while (0)
#define RS_ADDMUL(a, n, src, c, fl) \\
    do { mp_limb_t __cy = mpn_addmul_1(a, src, n, c); \\
         (a)[n] = (fl) + __cy; (fl) += ((a)[n] < __cy); } while (0)
/* first term of a block, negative sign: a[0..n] = -c*src wrapped */
#define RS_MUL_NEG(a, n, src, c, fl) \\
    do { (a)[n] = mpn_mul_1(a, src, n, c); \\
         (fl) = -(mp_limb_t) mpn_neg(a, a, (n) + 1); } while (0)
/* subtraction from a nonnegative accumulator (fl becomes live) */
#define RS_SUBMUL_Z(a, n, src, c, fl) \\
    do { mp_limb_t __cy = mpn_submul_1(a, src, n, c); \\
         (a)[n] = -__cy; (fl) = -(mp_limb_t) (__cy != 0); } while (0)
/* Accumulating variants, used where the window is full (g == w, i.e.
   zbits = 16): the slot a[n] is not fresh -- it already holds the
   previous term's carry limb -- so the borrow/carry is folded into
   it rather than assigned over it. */
#define RS_SUBMUL_ACC(a, n, src, c, fl) \\
    do { mp_limb_t __cy = mpn_submul_1(a, src, n, c); \\
         (fl) -= ((a)[n] < __cy); (a)[n] -= __cy; } while (0)
#define RS_ADDMUL_ACC(a, n, src, c, fl) \\
    do { mp_limb_t __cy = mpn_addmul_1(a, src, n, c); \\
         (a)[n] += __cy; (fl) += ((a)[n] < __cy); } while (0)
#define RS_SUBMUL_Z_ACC(a, n, src, c, fl) \\
    do { mp_limb_t __cy = mpn_submul_1(a, src, n, c); \\
         (fl) = -(mp_limb_t) ((a)[n] < __cy); (a)[n] -= __cy; } while (0)
""")
    return "\n".join(out)


class SGen:
    def __init__(self, func, xn, m, dr, zbits, name):
        assert func in ("atan", "atanh", "sincos", "sinhcosh")
        self.func, self.xn, self.m, self.dr = func, xn, m, dr
        self.zbits, self.name = zbits, name
        self.inl_declared = False
        self.Lz = 2 * zbits // 64
        self.alt = func in ("atan", "sincos")
        Nodd = max(1, (64 * xn) // (2 * zbits))
        Nev = (64 * xn + zbits) // (2 * zbits)
        if func in ("atan", "atanh"):
            self.NA = min(Nodd, ATAN_NMAX)
            self.sers = [] if self.NA == 1 else \
                [dict(kind="odd", K=list(range(0 if self.alt else 1,
                                               self.NA)),
                      cs="atan_series_cs", out="res", guard=None)]
        else:
            self.NS = min(Nodd, SC_NMAX)
            self.NC = min(Nev, SC_NMAX + 1)
            self.sers = []
            if self.NS >= 2:
                self.sers.append(dict(kind="odd",
                    K=list(range(0 if self.alt else 1, self.NS)),
                    cs="sin_series_cs", out="ysin",
                    guard="ysin != NULL"))
            if self.NC >= 3:
                self.sers.append(dict(kind="even",
                    K=list(range(2, self.NC)),
                    cs="cos_series_cs", out="ycos",
                    guard="ycos != NULL"))
        self.Nmax = max((s["K"][-1] + 1 for s in self.sers), default=1)
        if self.sers:
            assert m % 2 == 0 and m >= 2


    def inl_scratch(self, decl):
        if not self.inl_declared:
            decl.append("    mp_limb_t il0, il1, il2, il3;")
            self.inl_declared = True

    def inl_sqrhigh(self, body, r, a, size):
        """r[0..size-1] ~ high limbs of (a,size)^2, error < 3 ulp"""
        if size == 1:
            body.append("    umul_ppmm(%s[0], il0, %s[0], %s[0]);"
                        % (r, a, a))
        else:
            body.append("    umul_ppmm(%s[1], %s[0], %s[1], %s[1]);"
                        % (r, r, a, a))
            body.append("    umul_ppmm(il1, il0, %s[1], %s[0]);" % (a, a))
            body.append("    add_ssaaaa(%s[1], %s[0], %s[1], %s[0], "
                        "UWORD(0), (il1 << 1) | (il0 >> (FLINT_BITS - 1)));"
                        % (r, r, r, r))

    def inl_mulhigh2(self, body, r, a, b):
        """r[0..1] ~ high 2 limbs of (a,2)*(b,2), error < 3 ulp"""
        body.append("    umul_ppmm(%s[1], %s[0], %s[1], %s[1]);"
                    % (r, r, a, b))
        body.append("    umul_ppmm(il1, il0, %s[1], %s[0]);" % (a, b))
        body.append("    add_ssaaaa(%s[1], %s[0], %s[1], %s[0], "
                    "UWORD(0), il1);" % (r, r, r, r))
        body.append("    umul_ppmm(il1, il0, %s[0], %s[1]);" % (a, b))
        body.append("    add_ssaaaa(%s[1], %s[0], %s[1], %s[0], "
                    "UWORD(0), il1);" % (r, r, r, r))

    def inl_mul_1(self, body, r, a, size, c):
        """exact (r,size+1) = (a,size) * c"""
        if size == 1:
            body.append("    umul_ppmm(%s[1], %s[0], %s[0], %s);"
                        % (r, r, a, c))
        else:
            body.append("    umul_ppmm(il1, %s[0], %s[0], %s);"
                        % (r, a, c))
            body.append("    umul_ppmm(%s[2], %s[1], %s[1], %s);"
                        % (r, r, a, c))
            body.append("    add_ssaaaa(%s[2], %s[1], %s[2], %s[1], "
                        "UWORD(0), il1);" % (r, r, r, r))

    def inl_addmul_1(self, body, r, a, size, c):
        """exact (r,size+1) += (a,size) * c, top slot occupied"""
        if size == 1:
            body.append("    umul_ppmm(il1, il0, %s[0], %s);" % (a, c))
            body.append("    add_ssaaaa(%s[1], %s[0], %s[1], %s[0], "
                        "il1, il0);" % (r, r, r, r))
        else:
            body.append("    umul_ppmm(il1, il0, %s[0], %s);" % (a, c))
            body.append("    umul_ppmm(il3, il2, %s[1], %s);" % (a, c))
            body.append("    add_sssaaaaaa(%s[2], %s[1], %s[0], "
                        "%s[2], %s[1], %s[0], il3, il2, il0);"
                        % (r, r, r, r, r, r))
            body.append("    add_ssaaaa(%s[2], %s[1], %s[2], %s[1], "
                        "UWORD(0), il1);" % (r, r, r, r))

    def inl_boundary3(self, body, r, a, b):
        """r[0..2] ~ high 3 limbs of (a,3) * (0, b[0], b[1]),
        error < 4 ulp (used for the zero-padded z^m boundary)"""
        body.append("    umul_ppmm(%s[2], %s[1], %s[2], %s[1]);"
                    % (r, r, a, b))
        body.append("    umul_ppmm(il1, %s[0], %s[1], %s[1]);"
                    % (r, a, b))
        body.append("    add_ssaaaa(%s[2], %s[1], %s[2], %s[1], "
                    "UWORD(0), il1);" % (r, r, r, r))
        body.append("    umul_ppmm(il1, il0, %s[2], %s[0]);" % (a, b))
        body.append("    add_sssaaaaaa(%s[2], %s[1], %s[0], "
                    "%s[2], %s[1], %s[0], UWORD(0), il1, il0);"
                    % (r, r, r, r, r, r))
        body.append("    umul_ppmm(il1, il0, %s[0], %s[1]);" % (a, b))
        body.append("    add_ssaaaa(%s[1], %s[0], %s[1], %s[0], "
                    "UWORD(0), il1);" % (r, r, r, r))
        body.append("    umul_ppmm(il1, il0, %s[1], %s[0]);" % (a, b))
        body.append("    add_ssaaaa(%s[1], %s[0], %s[1], %s[0], "
                    "UWORD(0), il1);" % (r, r, r, r))

    def zl(self, j):
        """zero limbs of z^j.  For zbits = 16 the shrink is disabled:
        the powers' 32 j leading zero bits sit at half-limb offsets,
        which the limb-aligned power computation and window layout do
        not support (a shift-normalized z would be needed); the
        carry-slot collisions of full windows (g == w) are handled by
        the additive emissions below."""
        if self.zbits == 16:
            return 0
        return (2 * self.zbits * j) // 64

    def gn(self, j):
        return self.xn - self.zl(j)

    def src(self, p, g):
        off = self.gn(p) - g
        assert off >= 0
        name = "z" if p == 1 else "z%d" % p
        return name if off == 0 else "%s + %d" % (name, off)

    def dname(self):
        return ("ATAN_SERIES_D" if self.func in ("atan", "atanh")
                else "SIN_COS_SERIES_D")

    def iname(self):
        return ("atan_series_dinv" if self.func in ("atan", "atanh")
                else "sin_cos_series_dinv")

    def ilen(self):
        return INVA_LEN if self.func in ("atan", "atanh") else INVF_LEN

    def power_closure(self):
        m = self.m
        need = set()
        for s in self.sers:
            for k in s["K"]:
                p = k - (k // m) * m
                if p >= 2:
                    need.add(p)
            if s["K"][-1] >= m:
                need.add(m)
        stack = list(need)
        while stack:
            j = stack.pop()
            if j <= 2:
                continue
            a = j // 2
            for d in ((a, j - a) if j % 2 else (a,)):
                if d >= 2 and d not in need:
                    need.add(d)
                    stack.append(d)
        return sorted(need)

    def emit_series(self, body, S, bufs, decl=None):
        """full splitting for one series; returns (curbuf, othbuf, occ)"""
        m, xn, alt = self.m, self.xn, self.alt
        K, cs = S["K"], S["cs"]
        bases = sorted(set((k // m) * m for k in K), reverse=True)
        cur, oth = bufs
        occ, dyn = None, False
        for b in bases:
            w = xn - self.zl(b)
            for p in range(min(m - 1, K[-1] - b), -1, -1):
                k = b + p
                if k not in K:
                    continue
                g = xn - self.zl(b + p)
                c = "%s[%d]" % (cs, k)
                neg = alt and (k % 2 == 1)
                if p == 0:
                    if occ is not None and occ == w:
                        # zbits = 16 only: a previous carry limb sits
                        # at the units slot; add instead of assign
                        # (the window value keeps D < 2^63 headroom).
                        # This is correct with a live sign extension
                        # too: in wrapped arithmetic the add commits
                        # c * 2^(64 w) to the two's-complement value,
                        # and the discarded carry out of limb w is
                        # exactly what cancels fl_s = -1.  The block
                        # value is nonnegative (the units coefficient
                        # dominates), so the sign extension dies here.
                        assert not neg      # even m => even units k
                        body.append("    %s[%d] += %s;" % (cur, w, c))
                        dyn = False
                        continue
                    for i in range((0 if occ is None else occ + 1), w):
                        body.append("    %s[%d] = %s;"
                                    % (cur, i,
                                       "fl_s" if dyn else "0"))
                    body.append("    %s[%d] = %s%s;"
                                % (cur, w, c, " + fl_s" if dyn else ""))
                    occ, dyn = w, False
                    continue
                if occ is None:
                    if neg:
                        body.append("    RS_MUL_NEG(%s, %d, %s, %s, fl_s);"
                                    % (cur, g, self.src(p, g), c))
                        dyn = True
                    elif self.zbits == 16 and g <= 2:
                        self.inl_scratch(decl)
                        self.inl_mul_1(body, cur, self.src(p, g), g, c)
                    else:
                        body.append("    %s[%d] = mpn_mul_1(%s, %s, %d, %s);"
                                    % (cur, g, cur, self.src(p, g), g, c))
                    occ = g
                    continue
                for i in range(occ + 1, g):
                    body.append("    %s[%d] = %s;"
                                % (cur, i, "fl_s" if dyn else "0"))
                acc = (occ == g)        # full window: top slot in use
                if neg and dyn:
                    body.append("    RS_SUBMUL%s(%s, %d, %s, %s, fl_s);"
                                % ("_ACC" if acc else "",
                                   cur, g, self.src(p, g), c))
                elif neg:
                    body.append("    RS_SUBMUL_Z%s(%s, %d, %s, %s, fl_s);"
                                % ("_ACC" if acc else "",
                                   cur, g, self.src(p, g), c))
                    dyn = True
                elif dyn:
                    body.append("    RS_ADDMUL%s(%s, %d, %s, %s, fl_s);"
                                % ("_ACC" if acc else "",
                                   cur, g, self.src(p, g), c))
                elif occ == g:
                    # zbits = 16 only: the previous term's carry limb
                    # occupies cur[g]
                    if g <= 2:
                        self.inl_scratch(decl)
                        self.inl_addmul_1(body, cur, self.src(p, g), g, c)
                    else:
                        body.append(
                            "    %s[%d] += mpn_addmul_1(%s, %s, %d, %s);"
                            % (cur, g, cur, self.src(p, g), g, c))
                else:
                    body.append("    %s[%d] = mpn_addmul_1(%s, %s, %d, %s);"
                                % (cur, g, cur, self.src(p, g), g, c))
                occ = g
            if b > 0:
                assert not dyn        # even m: units precede boundaries
                n = w + 1
                gn_m = self.gn(m)
                zsrc = (self.src(m, n) if n <= gn_m else "z%d - 1" % m)
                if n > gn_m:
                    # z^m padded with one zero limb below; with the
                    # zbits = 16 full windows this occurs at every
                    # boundary, not only the lowest one
                    assert n == gn_m + 1
                if self.zbits == 16 and n == 3 and n > gn_m:
                    self.inl_scratch(decl)
                    self.inl_boundary3(body, oth, cur, "z%d" % m)
                else:
                    body.append("    flint_mpn_mulhigh_n(%s, %s, %s, %d);"
                                % (oth, cur, zsrc, n))
                cur, oth = oth, cur
                occ = n - 1
        return cur, oth, occ

    def emit_div(self, body, cur, oth, dlen, sub1):
        """q = floor(acc/D): in place (divrem) or into oth (inverse)"""
        if self.dr == 0:
            if sub1:
                body.append("    mpn_sub_1(%s, %s, %d, 1);"
                            % (cur, cur, dlen))
            body.append("    mpn_divrem_1(%s, 0, %s, %d, %s);"
                        % (cur, cur, dlen, self.dname()))
            return cur, oth
        off = self.ilen() - dlen
        if self.zbits == 16 and dlen == 2:
            self.inl_scratch(self._decl)
            self.inl_mulhigh2(body, oth, cur,
                              "(%s%s)" % (self.iname(),
                                          " + %d" % off if off else ""))
            return oth, cur
        body.append("    flint_mpn_mulhigh_n(%s, %s, %s%s, %d);"
                    % (oth, cur, self.iname(),
                       " + %d" % off if off else "", dlen))
        return oth, cur

    def head_even(self, body, qacc, qlen):
        """ycos = 1 -+ z/2 (+ q); the q addition cannot carry past
        limb gn_z - 1 since z/2 +- q < z"""
        xn, gn_z = self.xn, self.gn(1)
        body.append("    mpn_rshift(th, z, %d, 1);" % gn_z)
        if self.alt:
            body.append("    cy = mpn_neg(ycos, th, %d);" % gn_z)
            for i in range(gn_z, xn):
                body.append("    ycos[%d] = -cy;" % i)
            body.append("    ycos[%d] = 1 - cy;" % xn)
        else:
            body.append("    flint_mpn_copyi(ycos, th, %d);" % gn_z)
            for i in range(gn_z, xn):
                body.append("    ycos[%d] = 0;" % i)
            body.append("    ycos[%d] = 1;" % xn)
        if qacc is not None:
            body.append("    cy = mpn_add_n(ycos, ycos, %s, %d);"
                        % (qacc, qlen))
            body.append("    mpn_add_1(ycos + %d, ycos + %d, %d, cy);"
                        % (qlen, qlen, gn_z - qlen))

    def generate(self):
        func, xn, zbits = self.func, self.xn, self.zbits
        odd_only = func in ("atan", "atanh")
        if not self.sers:
            return self.generate_tiny()
        powers = self.power_closure()
        has_boundary = any(s["K"][-1] >= self.m for s in self.sers)
        need_head = (not odd_only) and self.NC >= 2
        gn_z = self.gn(1)

        decl, body = [], []
        self._decl = decl
        zdelta = 0 if zbits == 64 else self.zl(1)
        decl.append("    mp_limb_t z[%d];" % (gn_z + zdelta))
        for j in powers:
            size = self.gn(j)
            if j == self.m and has_boundary:
                decl.append("    mp_limb_t z%dp[%d];" % (j, size + 1))
            else:
                decl.append("    mp_limb_t z%d[%d];" % (j, size))
        decl.append("    mp_limb_t sa[%d];" % (xn + 1))
        decl.append("    mp_limb_t ta[%d];" % (xn + 1))
        if need_head:
            decl.append("    mp_limb_t th[%d];" % gn_z)

        pad = has_boundary and (self.m in powers)
        if pad:
            body.append("    mp_ptr z%d = z%dp + 1;" % (self.m, self.m))
            body.append("    z%dp[0] = 0;" % self.m)
            body.append("")
        if zbits == 64:
            body.append("    flint_mpn_sqrhigh(z, x + 1, %d);" % gn_z)
        elif zbits == 16 and xn <= 2:
            self.inl_scratch(decl)
            self.inl_sqrhigh(body, "z", "x", xn)
        else:
            body.append("    flint_mpn_sqrhigh(z, x, %d);" % xn)
        for j in powers:
            n = self.gn(j)
            a = j // 2
            if j % 2 == 0:
                if zbits == 16 and n <= 2:
                    self.inl_scratch(decl)
                    self.inl_sqrhigh(body, "z%d" % j, self.src(a, n), n)
                else:
                    body.append("    flint_mpn_sqrhigh(z%d, %s, %d);"
                                % (j, self.src(a, n), n))
            elif zbits == 16 and n == 2:
                self.inl_scratch(decl)
                self.inl_mulhigh2(body, "z%d" % j, self.src(a, n),
                                  self.src(a + 1, n))
            else:
                body.append("    flint_mpn_mulhigh_n(z%d, %s, %s, %d);"
                            % (j, self.src(a, n), self.src(a + 1, n), n))
        body.append("")

        for S in self.sers:
            sb, ind = [], ""
            cur, oth, occ = self.emit_series(sb, S, ("sa", "ta"),
                                             decl=decl)
            out = S["out"]
            xsig = xn - 1 if zbits == 64 else xn
            if S["kind"] == "odd" and self.alt:
                # q has xn significant limbs; res < 2^-zbits, so at
                # zbits = 64 its top limb is zero: shorten by one limb
                cur, oth = self.emit_div(sb, cur, oth, xn + 1, sub1=True)
                if zbits == 64:
                    sb.append("    flint_mpn_mulhigh_n(%s, %s + 1, x, %d);"
                              % (out, cur, xn - 1))
                    sb.append("    %s[%d] = 0;" % (out, xn - 1))
                else:
                    sb.append("    flint_mpn_mulhigh_n(%s, %s, x, %d);"
                              % (out, cur, xn))
            elif S["kind"] == "odd":
                # q = (S-1)*2^(64 xn) < 2^(64(xn-Lz)-2.5): its top
                # buffer limb is zero and x*(S-1) has n2 significant
                # limbs; the mpn_add_1 both copies x's upper limbs and
                # absorbs the carry (it cannot escape: x + x*(S-1) < 1)
                cur, oth = self.emit_div(sb, cur, oth, occ + 1, sub1=False)
                n2 = occ + xsig - xn
                if zbits == 16 and n2 == 2:
                    self.inl_scratch(decl)
                    self.inl_mulhigh2(sb, out, cur, "x")
                elif zbits == 16 and n2 == 1:
                    self.inl_scratch(decl)
                    sb.append("    umul_ppmm(%s[0], il0, %s[0], x[0]);"
                              % (out, cur))
                else:
                    sb.append("    flint_mpn_mulhigh_n(%s, %s + %d, "
                              "x + %d, %d);"
                              % (out, cur, occ - n2, xsig - n2, n2))
                if xsig - n2 == 0:
                    # zbits = 16: q x covers all limbs; the carry
                    # cannot escape (x + x (S-1) < 2^-15)
                    if n2 == 1:
                        sb.append("    %s[0] += x[0];" % out)
                    elif n2 == 2:
                        sb.append("    add_ssaaaa(%s[1], %s[0], %s[1], "
                                  "%s[0], x[1], x[0]);"
                                  % (out, out, out, out))
                    else:
                        sb.append("    mpn_add_n(%s, %s, x, %d);"
                                  % (out, out, n2))
                elif zbits == 64:
                    sb.append("    cy = mpn_add_n(%s, %s, x, %d);"
                              % (out, out, n2))
                    # x + x*(S-1) may cross 2^-64 when x is within an
                    # ulp of it: the carry-out is the top limb
                    sb.append("    %s[%d] = mpn_add_1(%s + %d, x + %d, "
                              "%d, cy);" % (out, xn - 1, out, n2, n2,
                                            xsig - n2))
                else:
                    sb.append("    cy = mpn_add_n(%s, %s, x, %d);"
                              % (out, out, n2))
                    sb.append("    mpn_add_1(%s + %d, x + %d, %d, cy);"
                              % (out, n2, n2, xsig - n2))
            else:
                cur, oth = self.emit_div(sb, cur, oth, occ + 1, sub1=False)
                self.head_even(sb, cur, occ)
            if out == "ysin":
                sb.append("    %s[%d] = 0;" % (out, xn))
            if S["guard"]:
                body.append("    if (%s)" % S["guard"])
                body.append("    {")
                body.extend("    " + ln for ln in sb)
                body.append("    }")
            else:
                body.extend(sb)
            body.append("")

        if not odd_only:
            extra = []
            if self.NS == 1:
                extra.append(("ysin != NULL",
                              ["    ysin[%d] = %s;"
                               % (i, "0" if (zbits == 64 and i == xn - 1)
                                  else "x[%d]" % i)
                               for i in range(xn)]
                              + ["    ysin[%d] = 0;" % xn]))
            if self.NC == 2 and not any(s["kind"] == "even"
                                        for s in self.sers):
                hb = []
                self.head_even(hb, None, 0)
                extra.append(("ycos != NULL", hb))
            elif self.NC == 1:
                extra.append(("ycos != NULL",
                              ["    ycos[%d] = 0;" % i for i in range(xn)]
                              + ["    ycos[%d] = 1;" % xn]))
            for guard, lines in extra:
                body.append("    if (%s)" % guard)
                body.append("    {")
                body.extend("    " + ln for ln in lines)
                body.append("    }")
                body.append("")

        bt = "\n".join(body)
        if "fl_s" in bt:
            decl.append("    mp_limb_t fl_s;")
            # With full windows (zbits = 16) the running sign
            # extension is consumed by the wrapped units addition
            # rather than read back, so its final value is dead:
            # say so, or the compiler rightly warns.
            if ("= fl_s" not in bt) and ("+ fl_s" not in bt):
                body.append("    (void) fl_s;")
                bt = "\n".join(body)
        if "cy" in bt:
            decl.append("    mp_limb_t cy;")
        while body and body[-1] == "":
            body.pop()
        eb = 12 if odd_only else 10
        return self.frame(decl, body, eb), eb

    def generate_tiny(self):
        func, xn, zbits = self.func, self.xn, self.zbits
        odd_only = func in ("atan", "atanh")
        decl, body = [], []
        if odd_only:
            for i in range(xn):
                body.append("    res[%d] = %s;"
                            % (i, "0" if (zbits == 64 and i == xn - 1)
                               else "x[%d]" % i))
            return self.frame([], body, 2), 2
        gn_z = self.gn(1)
        need_head = (self.NC >= 2)
        if need_head:
            zdelta = 0 if zbits == 64 else self.zl(1)
            decl.append("    mp_limb_t z[%d];" % (gn_z + zdelta))
            decl.append("    mp_limb_t th[%d];" % gn_z)
            if zbits == 64:
                body.append("    flint_mpn_sqrhigh(z, x + 1, %d);" % gn_z)
            else:
                body.append("    flint_mpn_sqrhigh(z, x, %d);" % xn)
        body.append("    if (ysin != NULL)")
        body.append("    {")
        for i in range(xn):
            body.append("        ysin[%d] = %s;"
                        % (i, "0" if (zbits == 64 and i == xn - 1)
                           else "x[%d]" % i))
        body.append("        ysin[%d] = 0;" % xn)
        body.append("    }")
        body.append("    if (ycos != NULL)")
        body.append("    {")
        hb = []
        if need_head:
            self.head_even(hb, None, 0)
        else:
            hb = ["    ycos[%d] = 0;" % i for i in range(xn)]
            hb.append("    ycos[%d] = 1;" % xn)
        body.extend("    " + ln for ln in hb)
        body.append("    }")
        bt = "\n".join(body)
        if "cy" in bt:
            decl.append("    mp_limb_t cy;")
        return self.frame(decl, body, 3 if need_head else 1), \
            (3 if need_head else 1)

    def frame(self, decl, body, eb):
        if not any("cy" in ln for ln in body):
            decl = [d for d in decl if d.strip() != "mp_limb_t cy;"]
        if any(d.strip() == "mp_limb_t il0, il1, il2, il3;"
               for d in decl):
            used = [v for v in ("il0", "il1", "il2", "il3")
                    if any(v in ln for ln in body)]
            decl = [d for d in decl
                    if d.strip() != "mp_limb_t il0, il1, il2, il3;"]
            if used:
                decl.append("    mp_limb_t %s;" % ", ".join(used))
        func, xn, zbits = self.func, self.xn, self.zbits
        pm = "-" if self.alt else "+"
        code = []
        if func in ("atan", "atanh"):
            code.append("/* res[0..%d] ~= sum_{k<%d} (%s1)^k x^(2k+1)/(2k+1);"
                        % (xn - 1, self.NA, pm))
            sig = "%s(mp_ptr res, mp_srcptr x)" % self.name
        else:
            o = ("ysin", "ycos") if func == "sincos" else ("ysinh", "ycosh")
            code.append("/* %s[0..%d] ~= sum_{k<%d} (%s1)^k x^(2k+1)/(2k+1)!,"
                        % (o[0], xn, self.NS, pm))
            code.append("   %s[0..%d] ~= sum_{k<%d} (%s1)^k x^(2k)/(2k)!;"
                        % (o[1], xn, self.NC, pm))
            code.append("   either output may be NULL;")
            sig = "%s(mp_ptr ysin, mp_ptr ycos, mp_srcptr x)" % self.name
        code.append("   x: %d frac limbs, 0 <= x < 2^-%d%s;"
                    % (xn, zbits,
                       " (top limb never read)" if zbits == 64 else ""))
        code.append("   %s error <= %d ulp (2^-%d) */"
                    % ("two-sided" if self.alt else "one-sided",
                       eb, 64 * xn))
        code.append("void")
        code.append(sig)
        code.append("{")
        code.extend(decl)
        if decl:
            code.append("")
        code.extend(body)
        code.append("}")
        return "\n".join(code)


def xnrange(func, zbits):
    if zbits == 16:
        # The shrink is disabled here (see SGen.zl), so every window
        # is full; the carry-slot collisions this creates are handled
        # by the accumulating RS_*_ACC emissions, which lets the
        # alternating families in as well.
        #
        # The reach is limited by the SHARED DENOMINATOR, not by the
        # emission: the series are clamped to ATAN_NMAX resp. SC_NMAX
        # terms, while x < 2^-16 needs about 2 xn of them.  The first
        # dropped term bounds the achievable precision:
        #   sin/cos: x^21/21! < 2^-401, enough for 64 xn <= 384
        #   atan:    x^35/35  < 2^-565, enough well past xn = 7
        # so the sin/cos families stop at xn = 6 (xn = 7 would carry
        # an error of some 2^46 ulp) and the atan families at xn = 7.
        return range(1, 8) if func in ("atan", "atanh") else range(1, 7)
    if func in ("atan", "atanh"):
        return range(1, (35 if zbits == 64 else 17) + 1)
    return range(1, (22 if zbits == 64 else 10) + 1)


def cfgrange(func, xn, zbits):
    g = SGen(func, xn, 2, 0, zbits, "probe")
    if not g.sers:
        return [(0, 0)]
    N = g.Nmax
    ms = [m for m in (2, 4, 6, 8) if m <= max(2, N - 1)] or [2]
    return [(m, dr) for m in ms for dr in (0, 1)]


def header():
    return """/* Generated by gen_atan_sincos.py -- do not edit.

   Hardcoded fixed-point evaluators for the atan/atanh and
   sin/cos/sinh/cosh Taylor series, 0 <= x < 2^-zbits. */

#include "longlong.h"
#include "mpn_extras.h"

"""


FUNCS = ("atan", "atanh", "sincos", "sinhcosh")
BASENAMES = {"atan": "mpn_atan_series", "atanh": "mpn_atanh_series",
             "sincos": "mpn_sin_cos_series",
             "sinhcosh": "mpn_sinh_cosh_series"}


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode.endswith("16"):
        zbits = 16
    elif mode.endswith("32"):
        zbits = 32
    else:
        zbits = 64
    mode = mode.replace("32", "").replace("16", "")
    sfx = "" if zbits == 64 else str(zbits)

    out = [header(), emit_tables()]
    out.append("typedef void (*odd_series_fn)(mp_ptr, mp_srcptr);")
    out.append("typedef void (*odd_even_series_fn)(mp_ptr, mp_ptr, mp_srcptr);")
    out.append("")

    if mode == "all":
        ent = {f: [] for f in FUNCS}
        for func in FUNCS:
            base = BASENAMES[func] + sfx
            for xn in xnrange(func, zbits):
                for (m, dr) in cfgrange(func, xn, zbits):
                    name = "%s_%d%s" % (base, xn,
                                        "_rs%d_d%d" % (m, dr) if m else "")
                    g = SGen(func, xn, max(m, 2), dr, zbits, name)
                    code, eb = g.generate()
                    out.append(code)
                    out.append("")
                    ent[func].append((name, xn, m, dr, eb))
        out.append("typedef struct { odd_series_fn f; int xn; int zbits;")
        out.append("    int m; int dr; int alternating; int errbound; }")
        out.append("odd_series_entry;")
        out.append("typedef struct { odd_even_series_fn f; int xn; int zbits;")
        out.append("    int m; int dr; int alternating; int errbound; }")
        out.append("odd_even_series_entry;")
        out.append("")
        out.append("odd_series_entry odd_series%s_all[] = {" % sfx)
        for func in ("atan", "atanh"):
            for (nm, xn, m, dr, eb) in ent[func]:
                out.append("    { %s, %d, %d, %d, %d, %d, %d },"
                           % (nm, xn, zbits, m, dr,
                              1 if func == "atan" else 0, eb))
        out.append("};")
        out.append("int odd_series%s_all_count = %d;"
                   % (sfx, len(ent["atan"]) + len(ent["atanh"])))
        out.append("")
        out.append("odd_even_series_entry odd_even_series%s_all[] = {" % sfx)
        for func in ("sincos", "sinhcosh"):
            for (nm, xn, m, dr, eb) in ent[func]:
                out.append("    { %s, %d, %d, %d, %d, %d, %d },"
                           % (nm, xn, zbits, m, dr,
                              1 if func == "sincos" else 0, eb))
        out.append("};")
        out.append("int odd_even_series%s_all_count = %d;"
                   % (sfx, len(ent["sincos"]) + len(ent["sinhcosh"])))
    elif mode == "best":
        with open(sys.argv[2]) as f:
            raw = json.load(f)
        for func in FUNCS:
            base = BASENAMES[func] + sfx
            xmax = max(xnrange(func, zbits))
            ebs = {}
            for xn in xnrange(func, zbits):
                m, dr = raw[func].get(str(xn), [2, 0])
                g = SGen(func, xn, max(int(m), 2), int(dr), zbits,
                         "%s_%d" % (base, xn))
                code, eb = g.generate()
                out.append(code)
                out.append("")
                ebs[xn] = eb
            tname = ("odd_series_fn" if func in ("atan", "atanh")
                     else "odd_even_series_fn")
            out.append("const %s %s_tab[%d] = {" % (tname, base, xmax + 1))
            for xn in range(xmax + 1):
                out.append("    %s," % ("NULL" if xn < 1
                                        else "%s_%d" % (base, xn)))
            out.append("};")
            out.append("const int %s_err_tab[%d] = {" % (base, xmax + 1))
            for xn in range(xmax + 1):
                out.append("    %d," % (0 if xn < 1 else ebs[xn]))
            out.append("};")
            out.append("")
    else:
        raise SystemExit("mode must be all[32] or best[32]")

    print("\n".join(out))




import re

# ===========================================================================
# Driver for the fixed module: regenerate src/fixed/trig_rs_hard.inc.
#
# The rectangular splitting parameters (m, use_divrem) below were
# selected by benchmarking every candidate on x86-64 (Skylake); to
# retune, generate mode "all", time the variants, and update the
# tables.
# ===========================================================================

WINNERS_64 = {'atanh': {1: (0, 0), 2: (0, 0), 3: (0, 0), 4: (2, 1), 5: (2, 1), 6: (2, 1), 7: (2, 1), 8: (2, 1), 9: (2, 1), 10: (4, 1), 11: (4, 1), 12: (4, 0), 13: (4, 0), 14: (4, 0), 15: (6, 0), 16: (4, 0), 17: (6, 0), 18: (6, 0), 19: (6, 0), 20: (6, 0), 21: (4, 0), 22: (4, 0), 23: (4, 0), 24: (4, 0), 25: (4, 0), 26: (4, 0), 27: (4, 0), 28: (4, 0), 29: (4, 0), 30: (4, 0), 31: (4, 0), 32: (4, 0), 33: (4, 0), 34: (4, 0), 35: (4, 0)}, 'atan': {1: (0, 0), 2: (0, 0), 3: (0, 0), 4: (2, 1), 5: (2, 1), 6: (2, 1), 7: (2, 1), 8: (2, 1), 9: (2, 0), 10: (4, 0), 11: (4, 0), 12: (4, 0), 13: (4, 0), 14: (4, 0), 15: (6, 0), 16: (4, 0), 17: (2, 0), 18: (6, 0), 19: (2, 0), 20: (6, 0), 21: (4, 0), 22: (4, 0), 23: (4, 0), 24: (4, 0), 25: (4, 0), 26: (4, 0), 27: (4, 0), 28: (6, 0), 29: (4, 0), 30: (4, 0), 31: (4, 0), 32: (4, 0), 33: (4, 0), 34: (4, 0), 35: (4, 0)}, 'sinhcosh': {1: (0, 0), 2: (0, 0), 3: (0, 0), 4: (2, 1), 5: (2, 1), 6: (2, 1), 7: (2, 1), 8: (2, 1), 9: (4, 1), 10: (4, 1), 11: (4, 1), 12: (4, 1), 13: (6, 0), 14: (6, 0), 15: (6, 0), 16: (6, 0), 17: (8, 0), 18: (8, 0), 19: (8, 0), 20: (8, 0), 21: (8, 0), 22: (8, 0)}, 'sincos': {1: (0, 0), 2: (0, 0), 3: (0, 0), 4: (2, 1), 5: (2, 1), 6: (2, 1), 7: (2, 1), 8: (2, 1), 9: (4, 1), 10: (4, 1), 11: (4, 0), 12: (4, 0), 13: (6, 0), 14: (6, 0), 15: (6, 0), 16: (6, 0), 17: (6, 0), 18: (8, 0), 19: (8, 0), 20: (6, 0), 21: (8, 0), 22: (8, 0)}}

WINNERS_16 = {'atanh': {1: (2, 1), 2: (2, 1), 3: (2, 1), 4: (4, 1),
              5: (6, 1), 6: (6, 1), 7: (6, 1)},
              'atan': {1: (2, 1), 2: (2, 1), 3: (2, 1), 4: (4, 1),
              5: (4, 1), 6: (4, 1), 7: (4, 1)},
              'sincos': {1: (2, 1), 2: (2, 1), 3: (2, 1), 4: (4, 1),
              5: (6, 1), 6: (6, 1)},
              'sinhcosh': {1: (2, 1), 2: (2, 1), 3: (2, 1), 4: (4, 1),
              5: (6, 1), 6: (6, 1)}}

WINNERS_32 = {'atanh': {1: (0, 0), 2: (2, 1), 3: (2, 1), 4: (2, 1), 5: (4, 1), 6: (4, 1), 7: (6, 1), 8: (4, 1), 9: (4, 1), 10: (4, 0), 11: (6, 0), 12: (4, 0), 13: (6, 0), 14: (4, 0), 15: (4, 0), 16: (4, 0), 17: (4, 0)}, 'atan': {1: (0, 0), 2: (2, 1), 3: (2, 1), 4: (2, 1), 5: (4, 1), 6: (4, 1), 7: (6, 1), 8: (4, 1), 9: (4, 1), 10: (6, 0), 11: (4, 0), 12: (4, 0), 13: (4, 0), 14: (4, 0), 15: (4, 0), 16: (4, 0), 17: (4, 0)}, 'sinhcosh': {1: (0, 0), 2: (2, 1), 3: (2, 1), 4: (2, 1), 5: (4, 1), 6: (4, 1), 7: (6, 1), 8: (6, 1), 9: (8, 1), 10: (8, 1)}, 'sincos': {1: (0, 0), 2: (2, 1), 3: (2, 1), 4: (2, 1), 5: (4, 1), 6: (4, 1), 7: (6, 1), 8: (6, 1), 9: (8, 1), 10: (8, 0)}}


def generate_best(zbits, best, funcs=None):
    sfx = "" if zbits == 64 else ("32" if zbits == 32 else "16")
    out = [header(), emit_tables()]
    out.append("typedef void (*odd_series_fn)(mp_ptr, mp_srcptr);")
    out.append("typedef void (*odd_even_series_fn)(mp_ptr, mp_ptr, mp_srcptr);")
    out.append("")
    for func in (funcs if funcs is not None else FUNCS):
        base = BASENAMES[func] + sfx
        xmax = max(xnrange(func, zbits))
        ebs = {}
        for xn in xnrange(func, zbits):
            m, dr = best[func].get(xn, (2, 0))
            g = SGen(func, xn, max(int(m), 2), int(dr), zbits,
                     "%s_%d" % (base, xn))
            code, eb = g.generate()
            out.append(code)
            out.append("")
            ebs[xn] = eb
        tname = ("odd_series_fn" if func in ("atan", "atanh")
                 else "odd_even_series_fn")
        out.append("const %s %s_tab[%d] = {" % (tname, base, xmax + 1))
        for xn in range(xmax + 1):
            out.append("    %s," % ("NULL" if xn < 1
                                    else "%s_%d" % (base, xn)))
        out.append("};")
        out.append("const int %s_err_tab[%d] = {" % (base, xmax + 1))
        for xn in range(xmax + 1):
            out.append("    %d," % (0 if xn < 1 else ebs[xn]))
        out.append("};")
        out.append("")
    return "\n".join(out)


def _convert_types(s):
    s = s.replace('mp_srcptr', 'nn_srcptr').replace('mp_ptr', 'nn_ptr')
    s = s.replace('mp_limb_t', 'ulong')
    return s


def _drop_err_tabs(s):
    s = re.sub(r'(static )?const int \w+_err_tab\[\d+\] = \{[^}]*\};\n*', '', s)
    return s


def _statify(s):
    s = re.sub(r'\nvoid\nmpn_(atan|atanh|sin_cos|sinh_cosh)_series',
               r'\nstatic void\nmpn_\1_series', s)
    s = s.replace('\nconst odd_series_fn ', '\nstatic const odd_series_fn ')
    s = s.replace('\nconst odd_even_series_fn ',
                  '\nstatic const odd_even_series_fn ')
    s = s.replace('typedef void (*odd_series_fn)(mp_ptr, mp_srcptr);', '')
    s = s.replace('typedef void (*odd_even_series_fn)(mp_ptr, mp_ptr, mp_srcptr);', '')
    return s


if __name__ == "__main__":
    # The 2^-16 window family died with the r >= 32 contract of the
    # public bitwise functions (the fully specialized per-size paths
    # in *_opt_<n>.c serve everything below that); only the 2^-64 and
    # 2^-32 families remain.
    f64 = generate_best(64, WINNERS_64)
    f32 = generate_best(32, WINNERS_32)

    f64 = f64[f64.index('/* atan/atanh: c_k'):]
    m = re.search(r'\nvoid\nmpn_\w+_series32_1\(', f32)
    f32 = f32[m.start():]

    out = ("/* Generated by dev/gen_fixed_trig_rs_hard.py -- private static\n"
           "   Taylor series routines for the fixed module -- do not edit"
           " by hand. */\n\n"
           "typedef void (*odd_series_fn)(nn_ptr, nn_srcptr);\n"
           "typedef void (*odd_even_series_fn)(nn_ptr, nn_ptr, nn_srcptr);\n\n"
           + _statify(f64) + "\n" + _statify(f32))
    out = _drop_err_tabs(out)
    for a, b in [('mpn_atanh_series16_', '_fixed_atanh_rs16_'),
                 ('mpn_atan_series16_', '_fixed_atan_rs16_'),
                 ('mpn_sin_cos_series16_', '_fixed_sin_cos_rs16_'),
                 ('mpn_atan_series32_', '_fixed_atan_rs32_'),
                 ('mpn_atanh_series32_', '_fixed_atanh_rs32_'),
                 ('mpn_sin_cos_series32_', '_fixed_sin_cos_rs32_'),
                 ('mpn_sinh_cosh_series32_', '_fixed_sinh_cosh_rs32_'),
                 ('mpn_atan_series_', '_fixed_atan_rs_'),
                 ('mpn_atanh_series_', '_fixed_atanh_rs_'),
                 ('mpn_sin_cos_series_', '_fixed_sin_cos_rs_'),
                 ('mpn_sinh_cosh_series_', '_fixed_sinh_cosh_rs_')]:
        out = out.replace(a, b)
    out = _convert_types(out)
    print(out.rstrip("\n"))
