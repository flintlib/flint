#!/usr/bin/env python3
"""Regenerate src/fixed/exp_rs_hard.inc (run from the top-level
FLINT directory: python3 dev/gen_fixed_exp_rs_hard.py > src/fixed/exp_rs_hard.inc)."""

#!/usr/bin/env python3
"""
Generator for hardcoded fixed-point exp Taylor evaluators, generalizing
mpn_exp_series_10_rs4 to both x < 2^-64 (zbits = 64) and x < 2^-32
(zbits = 32).

A generated function

    void mpn_exp_series[32]_{N}[_rs{m}](mp_ptr res, mp_srcptr x)

computes res (xn+1 limbs, xn fractional + 1 integer) ~= sum_{k<N} x^k/k!
for fixed-point x with xn fractional limbs, 0 <= x < 2^-zbits, where
xn = ceil(N*zbits/64) (so xn = N for zbits = 64, xn = ceil(N/2) for
zbits = 32; for zbits = 64 the top limb of x is never read).  Terms
1, x, x^2/2 are added exactly; terms k = 3..N-1 are evaluated by
rectangular splitting with block length m over the common denominator
D = 20! (coefficients 20!/k! fit one limb, hence N <= 21), with a
single final division by D as a mulhigh by a truncated inverse or as
mpn_divrem_1.

Everything is limb-granular via lz(j) = floor(j*zbits/64) leading zero
limbs of x^j:

  * x^j is stored grid-aligned (bottom limb at weight 2^(-64*xn)) with
    gn(j) = xn - lz(j) significant limbs; products use factor pairs
    (a, b) with lz(a) + lz(b) = lz(j) so the mulhigh length is exactly
    n = xn - lz(a) - lz(b) (for zbits = 32: odd j from (j//2, j//2+1),
    j = 0 mod 4 by squaring, j = 2 mod 4 from (j/2-1, j/2+1); x^2 is
    the one delta = 1 case, computed as a full sqrhigh of x with one
    extra always-zero top limb).
  * Inside a block with base b (window w = xn - lz(b)), term k = b + p
    reads the top g = w - lz(p) limbs of x^p.  For zbits = 64, g grows
    by one limb every term and each op's carry limb is stored fresh;
    for zbits = 32, g grows every other term and an op whose carry
    limb is already occupied adds into it instead -- at most two
    carries (< 2 * 20! < 2^64) meet in one limb, so this cannot
    overflow further.
  * Block bases must satisfy lz(b) - lz(m) = lz(b - m), which holds
    for all bases exactly when m is even (or zbits = 64); for
    zbits = 32 only even m is generated.
  * Boundaries multiply the alen occupied limbs by the top alen limbs
    of x^m; only the last boundary (base m -> 0) needs x^m
    zero-extended one limb below (padded buffer).
  * The RS sum is a pure fraction < c_3 x^3 * 1.1 in D-scaled grid
    units, so it occupies xn - 2 (zbits = 64) resp. xn (zbits = 32)
    limbs with headroom, and reading the top dlen limbs of the
    precomputed floor(2^(64*20)/20!) keeps the inverse-truncation
    error below 2^-5 ulp regardless of N.

Error bound (ulp = 2^(-64*xn), one-sided: computed <= exact):
  division 2 + 1/32; x^2/2 term <= 2.5; power errors weighted by
  1/k! <= 0.5; boundary x^m read errors <= 6.6/m!; window truncations
  are all damped by the final division by 20! (2^-61).  Total <= 8 ulp
  for m >= 3 and <= 11 ulp for m = 2.

Usage:
  python3 gen_exp_series.py all[32]  > file.c        # every (N, m, div)
  python3 gen_exp_series.py best[32] M_MAP_FILE > file.c
"""

import sys
import json

FAC20 = 1
for i in range(1, 21):
    FAC20 *= i

QLIMBS = 20
QINV = (2 ** (64 * QLIMBS)) // FAC20
assert QINV < 2 ** (64 * QLIMBS - 60)
assert QINV >= 2 ** (64 * (QLIMBS - 1))

NMAX = 21


def limbs(v, n):
    return [(v >> (64 * i)) & ((1 << 64) - 1) for i in range(n)]


def emit_tables():
    out = []
    out.append("/* c_k = 20!/k!, equal to factorial_tab_numer[0..20] in arb */")
    out.append("static const mp_limb_t exp_series_cs[21] = {")
    for k in range(21):
        c = FAC20
        for i in range(1, k + 1):
            c //= i
        out.append("    UWORD(%d)," % c)
    out.append("};")
    out.append("")
    out.append("/* floor(2^(64*%d) / 20!), %d significant limbs */" % (QLIMBS, QLIMBS))
    out.append("static const mp_limb_t exp_series_fac20_inv[%d] = {" % QLIMBS)
    for l in limbs(QINV, QLIMBS):
        out.append("    UWORD(%d)," % l)
    out.append("};")
    out.append("")
    return "\n".join(out)


class Gen:
    def __init__(self, N, m, zbits, use_divrem=False, name=None,
                 xn=None, noshrink=False):
        # zbits need not be a divisor of 64: the leading-zero shrink
        # lz(j) = floor(j zbits / 64) is well defined for any zbits,
        # and xn (the output precision) may be given explicitly rather
        # than derived from N, which is what lets the series spend the
        # 1/N! gain on FEWER TERMS instead of on slack: N is then the
        # smallest number of terms with N zbits + log2(N!) >= 64 xn.
        assert 4 <= zbits <= 64
        assert 1 <= N <= NMAX
        self.N = N
        self.m = m
        self.zbits = zbits
        self.noshrink = noshrink
        self.m_split = (m != 0 and N >= 4
                        and (zbits == 64 or m % (64 // zbits) == 0))
        self.use_divrem = use_divrem
        self.xn = xn if xn is not None else -(-N * zbits // 64)
        base = ("mpn_exp_series" if zbits == 64
                else "mpn_exp_series%d" % zbits)
        self.name = name or "%s_%d_rs%d%s" % (base, N, m,
                                              "_dr" if use_divrem else "")
        if N >= 4 and m:
            assert 2 <= m <= N
            # Block bases need lz(b) - lz(m) = lz(b - m).  With m = N
            # there is a single block (base 0) and no boundary, so the
            # constraint is vacuous -- which is what makes arbitrary
            # zbits usable at all.  Otherwise it forces zbits | 64 and
            # m a multiple of 64 / zbits.
            if m < N and not self.noshrink and (64 % zbits != 0
                          or m % (64 // zbits) != 0):
                raise ValueError("m = %d illegal at zbits = %d"
                                 % (m, zbits))

    def lz(self, j):
        # With the shrink on, lz(j) = floor(j zbits / 64) is the number
        # of leading zero limbs of x^j.  It satisfies the block-base
        # identity lz(b) - lz(m) = lz(b - m) only when 64 % zbits == 0,
        # which is exactly why the tuned families sit at zbits = 16, 32
        # and 64: at any other zbits rectangular splitting is illegal
        # and the series must be a single block.  Turning the shrink OFF
        # gives up the leading-zero savings but makes lz identically 0,
        # so the identity is trivial and ANY m is legal -- the tradeoff
        # worth measuring for bit-granular r.
        if self.noshrink:
            return 0
        return (j * self.zbits) // 64

    def gn(self, j):
        return self.xn - self.lz(j)

    def xsig(self):
        # significant limbs of x itself (= gn(1))
        return self.gn(1)

    def src(self, p, g):
        """top g limbs of x^p's grid part"""
        if p == 1:
            off = self.xsig() - g
            name = "x"
        else:
            off = self.gn(p) - g
            name = "x%d" % p
        assert off >= 0
        return name if off == 0 else "%s + %d" % (name, off)

    def pair(self, j):
        """factor pair (a, b), a <= b, a + b = j, minimizing the
           alignment delta lz(j) - lz(a) - lz(b) (0 except x^2 at
           zbits = 32, where delta = 1)"""
        best = None
        for a in range(j // 2, 0, -1):
            b = j - a
            d = self.lz(j) - self.lz(a) - self.lz(b)
            assert d >= 0
            key = (d, 0 if a == b else 1)
            if best is None or key < best[0]:
                best = (key, a, b)
        return best[1], best[2], best[0][0]

    def power_closure(self, term_powers, has_boundary):
        need = set(p for p in term_powers if p >= 2)
        if has_boundary:
            need.add(self.m)
        stack = list(need)
        while stack:
            j = stack.pop()
            if j <= 2:
                continue
            a, b, _ = self.pair(j)
            for d in (a, b):
                if d >= 2 and d not in need:
                    need.add(d)
                    stack.append(d)
        return sorted(need)

    def generate(self):
        N, m, xn, zbits = self.N, self.m, self.xn, self.zbits
        body, decl = [], []
        errbound = 11 if (N >= 5 and m == 2) else 8

        if N <= 3:
            return self.generate_tiny(errbound)

        bases = sorted(set((k // m) * m for k in range(3, N)), reverse=True)
        has_boundary = len(bases) > 1 or (bases and bases[0] != 0)
        term_powers = set(k - (k // m) * m for k in range(3, N))
        powers = self.power_closure(term_powers, has_boundary)

        # buffers
        for j in powers:
            _, _, delta = self.pair(j)
            size = self.gn(j) + delta
            if j == m and has_boundary:
                decl.append("    mp_limb_t x%dp[%d];" % (j, size + 1))
            else:
                decl.append("    mp_limb_t x%d[%d];" % (j, size))
        decl.append("    mp_limb_t s[%d];" % (xn + 1))
        decl.append("    mp_limb_t t[%d];" % (xn + 1))

        if has_boundary:
            body.append("    mp_ptr x%d = x%dp + 1;" % (m, m))
            body.append("    x%dp[0] = 0;" % m)
            body.append("")

        # powers
        for j in powers:
            a, b, delta = self.pair(j)
            n = xn - self.lz(a) - self.lz(b)
            assert n == self.gn(j) + delta
            assert n <= self.gn(a) and n <= self.gn(b)
            if a == b:
                body.append("    flint_mpn_sqrhigh(x%d, %s, %d);"
                            % (j, self.src(a, n), n))
            else:
                body.append("    flint_mpn_mulhigh_n(x%d, %s, %s, %d);"
                            % (j, self.src(a, n), self.src(b, n), n))
        body.append("")

        # rectangular splitting, descending k; alen = occupied acc limbs
        acc, other = "t", "s"
        first_op = True
        alen = 0
        for b in bases:
            w = xn - self.lz(b)
            khi = min(N - 1, b + m - 1)
            klo = max(3, b)
            for k in range(khi, klo - 1, -1):
                p = k - b
                g = w - self.lz(p)
                c = "exp_series_cs[%d]" % k
                if first_op:
                    if p == 0:
                        for i in range(w):
                            body.append("    %s[%d] = 0;" % (acc, i))
                        body.append("    %s[%d] = %s;" % (acc, w, c))
                        alen = w + 1
                    elif g == 1:
                        sname = "x" if p == 1 else "x%d" % p
                        stop = (self.xsig() - 1) if p == 1 else (self.gn(p) - 1)
                        body.append("    umul_ppmm(%s[1], %s[0], %s[%d], %s);"
                                    % (acc, acc, sname, stop, c))
                        alen = 2
                    else:
                        body.append("    %s[%d] = mpn_mul_1(%s, %s, %d, %s);"
                                    % (acc, g, acc, self.src(p, g), g, c))
                        alen = g + 1
                    first_op = False
                elif p == 0:
                    assert alen in (w, w + 1)
                    op = "+=" if alen == w + 1 else "="
                    body.append("    %s[%d] %s %s;" % (acc, w, op, c))
                    alen = w + 1
                else:
                    assert alen in (g, g + 1)
                    op = "+=" if alen == g + 1 else "="
                    body.append("    %s[%d] %s mpn_addmul_1(%s, %s, %d, %s);"
                                % (acc, g, op, acc, self.src(p, g), g, c))
                    alen = g + 1
            if b > 0:
                n = alen
                gn_m = self.gn(m)
                if n <= gn_m:
                    xm_src = self.src(m, n)
                else:
                    assert n == gn_m + 1 and b == m
                    xm_src = "x%d - 1" % m
                body.append("    flint_mpn_mulhigh_n(%s, %s, %s, %d);"
                            % (other, acc, xm_src, n))
                acc, other = other, acc
                alen = n
        body.append("")

        # final division by 20!
        dlen = alen
        if zbits == 64:
            assert dlen == xn - 2
        elif zbits == 32:
            assert dlen == xn
        if self.use_divrem:
            body.append("    mpn_divrem_1(res, 0, %s, %d, exp_series_cs[0]);"
                        % (acc, dlen))
        else:
            body.append("    flint_mpn_mulhigh_n(res, %s, exp_series_fac20_inv + %d, %d);"
                        % (acc, QLIMBS - dlen, dlen))
        body.append("")

        # + x^2/2, + x, + 1
        g2 = self.gn(2)
        body.append("    mpn_rshift(%s, x2, %d, 1);" % (other, g2))
        if dlen == g2:
            # carry limb of the x^2/2 add is fresh (zbits = 64)
            body.append("    res[%d] = mpn_add_n(res, res, %s, %d);"
                        % (g2, other, g2))
        else:
            # carry limb is the division output's top limb (zbits = 32);
            # no further carry: the total value is < 2^-31 < 1
            assert dlen == g2 + 1
            body.append("    res[%d] += mpn_add_n(res, res, %s, %d);"
                        % (g2, other, g2))
        xs = self.xsig()
        if xs + 1 <= xn:
            body.append("    res[%d] = mpn_add_n(res, res, x, %d);" % (xs, xs))
            body.append("    res[%d] = 1;" % xn)
        else:
            # x spans all xn limbs; the carry is provably zero
            # (total value < 2^-31), kept for robustness
            body.append("    res[%d] = 1 + mpn_add_n(res, res, x, %d);"
                        % (xn, xs))

        return self.frame(decl, body, errbound), errbound

    def generate_tiny(self, errbound):
        """N <= 3: no rectangular splitting part"""
        N, xn = self.N, self.xn
        body, decl = [], []
        if N == 1:
            for i in range(xn):
                body.append("    res[%d] = 0;" % i)
            body.append("    res[%d] = 1;" % xn)
        elif N == 2:
            xs = self.xsig()
            for i in range(xs):
                body.append("    res[%d] = x[%d];" % (i, i))
            for i in range(xs, xn):
                body.append("    res[%d] = 0;" % i)
            body.append("    res[%d] = 1;" % xn)
        else:  # N == 3: 1 + x + x^2/2
            decl.append("    mp_limb_t x2[%d];" % (self.gn(2) + self.pair(2)[2]))
            decl.append("    mp_limb_t t[%d];" % (self.gn(2) + 1))
            a, b, delta = self.pair(2)
            n = xn - self.lz(a) - self.lz(b)
            body.append("    flint_mpn_sqrhigh(x2, %s, %d);" % (self.src(1, n), n))
            g2 = self.gn(2)
            body.append("    mpn_rshift(t, x2, %d, 1);" % g2)
            for i in range(g2):
                body.append("    res[%d] = t[%d];" % (i, i))
            for i in range(g2, xn):
                body.append("    res[%d] = 0;" % i)
            xs = self.xsig()
            if xs + 1 <= xn:
                body.append("    res[%d] = mpn_add_n(res, res, x, %d);" % (xs, xs))
                body.append("    res[%d] = 1;" % xn)
            else:
                body.append("    res[%d] = 1 + mpn_add_n(res, res, x, %d);"
                            % (xn, xs))
        errbound = 3   # x^2/2: <= 2.5; N <= 2 exact
        return self.frame(decl, body, errbound), errbound

    def frame(self, decl, body, errbound):
        N, xn, zbits = self.N, self.xn, self.zbits
        code = []
        code.append("/* res[0..%d] ~= sum_{k<%d} x^k/k!; x: %d frac limbs,"
                    % (xn, N, xn))
        code.append("   0 <= x < 2^-%d; one-sided error <= %d ulp (2^-%d) */"
                    % (zbits, errbound, 64 * xn))
        code.append("void")
        code.append("%s(mp_ptr res, mp_srcptr x)" % self.name)
        code.append("{")
        code.extend(decl)
        if decl:
            code.append("")
        code.extend(body)
        code.append("}")
        return "\n".join(code)


# Hand-optimized zbits = 32 functions for N = 3, 4, 5, using inline limb
# operations.  Facts exploited: for xn = 2 the top limb x[1] < 2^32, so
# x[1]^2 is a plain single-limb multiply (no high word); x^2 * 2^(64*xn)
# fits few limbs with known headroom; floor(v/6) and floor(v/24) of a
# two-limb v with small high limb v1 are computed branch-free and
# one-sided as (v1 + ch)*floor(2^64/D) + floor(cl/D), where ch:cl =
# v0 + (2^64 mod D)*v1, underestimating by at most 1; and for N = 5,
# x^3/6 + x^4/24 = (4*x^3 + x^4)/24 needs only one division.
HAND32 = {}

HAND32[3] = ("""\
/* res[0..2] ~= 1 + x + x^2/2; x: 2 frac limbs, 0 <= x < 2^-32;
   one-sided error <= 3 ulp (2^-128) */
void
mpn_exp_series32_3%s(mp_ptr res, mp_srcptr x)
{
    mp_limb_t h, l, t;

    /* x^2 * 2^128 fits one limb; x[1] < 2^32 so x[1]^2 has no high word */
    umul_ppmm(h, l, x[1], x[0]);
    t = x[1] * x[1] + ((h << 1) | (l >> 63));
    add_ssaaaa(res[1], res[0], x[1], x[0], 0, t >> 1);
    res[2] = 1;
}""", 3)

HAND32[4] = ("""\
/* res[0..2] ~= 1 + x + x^2/2 + x^3/6; x: 2 frac limbs, 0 <= x < 2^-32;
   one-sided error <= 4 ulp (2^-128) */
void
mpn_exp_series32_4%s(mp_ptr res, mp_srcptr x)
{
    mp_limb_t h, l, s, h3, l3, t;

    s = x[1] * x[1];                    /* exact: x[1] < 2^32 */
    umul_ppmm(h3, l3, s, x[1]);         /* x^3 * 2^128 ~= h3 < 2^32 */
    umul_ppmm(h, l, x[1], x[0]);
    t = s + ((h << 1) | (l >> 63));     /* x^2 * 2^128 */
    t = (t >> 1) + h3 / 6;              /* x^2/2 + x^3/6 */
    add_ssaaaa(res[1], res[0], x[1], x[0], 0, t);
    res[2] = 1;
}""", 4)

HAND32[5] = ("""\
/* res[0..3] ~= 1 + x + x^2/2 + x^3/6 + x^4/24; x: 3 frac limbs,
   0 <= x < 2^-32; one-sided error <= 5 ulp (2^-192).
   All quantities are underestimates, so the result is one-sided. */
void
mpn_exp_series32_5%s(mp_ptr res, mp_srcptr x)
{
    mp_limb_t x0 = x[0], x1 = x[1], x2 = x[2];
    mp_limb_t ah, al, bh, bl, ch, cl, t1h, t1l, mh, ml, m2, c0, cA, cy;
    mp_limb_t b1, b0, d1, d0, e1, e0, f1, f0, v1, v0, w1, w0;
    mp_limb_t T1, T0, q1, q0, u, r;

    /* b1:b0 = floor(x^2 * 2^192): X^2 limbs [4:3], X = x2:x1:x0,
       dropping the sub-b0 part (<= 2 units) */
    umul_ppmm(ah, al, x2, x1);
    umul_ppmm(bh, bl, x2, x0);
    umul_ppmm(ch, cl, x1, x1);
    t1h = (ah << 1) | (al >> 63);       /* 2*x2*x1, high */
    t1l = al << 1;
    ml = cl + (bl << 1);                /* mid = x1^2 + 2*x2*x0 */
    c0 = ml < cl;
    mh = ch + ((bh << 1) | (bl >> 63));
    cA = mh < ch;
    mh += c0;
    m2 = cA | (mh < c0);
    b0 = t1l + mh;
    cy = b0 < t1l;
    b1 = x2 * x2 + t1h + m2 + cy;       /* x2 < 2^32: single multiply */

    /* v1:v0 = floor(x^3 * 2^192) = floor(X*b / 2^192) (drop <= 3 units);
       note e1 + f1 + carry can exceed 2^64, so fold with carried adds */
    umul_ppmm(d1, d0, x2, b1);
    umul_ppmm(e1, e0, x2, b0);
    umul_ppmm(f1, f0, x1, b1);
    add_ssaaaa(v1, v0, d1, d0, 0, e1);
    add_ssaaaa(v1, v0, v1, v0, 0, f1);
    add_ssaaaa(v1, v0, v1, v0, 0, (mp_limb_t)((e0 + f0) < e0));

    /* w1 = floor(x^4 * 2^192) = floor(b^2 / 2^192) (drop <= 3 units) */
    umul_ppmm(w1, w0, b1, b1);

    /* T = 4*x^3 + x^4, then floor(T/24) in one division:
       floor(T/24) >= (T1 + ch)*floor(2^64/24) + floor(cl/24),
       ch:cl = T0 + 16*T1 (2^64 mod 24 = 16), short by at most 1 */
    T0 = v0 << 2;
    T1 = (v1 << 2) | (v0 >> 62);
    T0 += w1;
    T1 += (T0 < w1);
    cl = T0 + 16 * T1;                  /* 16*T1 < 2^40 */
    ch = cl < T0;
    umul_ppmm(q1, q0, T1 + ch, UWORD(768614336404564650));
    u = cl / 24;
    q0 += u;
    q1 += (q0 < u);

    /* + x^2/2 */
    add_ssaaaa(q1, q0, q1, q0, b1 >> 1, (b0 >> 1) | (b1 << 63));

    /* res = 1 + x + q */
    r = x0 + q0;
    cy = r < x0;
    res[0] = r;
    r = x1 + q1;
    c0 = r < x1;
    r += cy;
    c0 += (r < cy);
    res[1] = r;
    res[2] = x2 + c0;
    res[3] = 1;
}""", 5)


def header():
    return """/* Generated by gen_exp_series.py -- do not edit.

   Hardcoded fixed-point evaluators for sum_{k<N} x^k/k!,
   xn = ceil(N*zbits/64) limbs, 0 <= x < 2^-zbits. */

#include "longlong.h"
#include "mpn_extras.h"

"""


def mrange(N, zbits):
    if N <= 3:
        return [0]                      # m irrelevant
    # at zbits = 16 the block size must be a multiple of 4, so allow
    # m == N as well (a single block: plain Horner over the terms),
    # otherwise N = 4 would admit no legal splitting at all
    hi = min(N if zbits == 16 else N - 1, 12)
    ms = range(2, hi + 1)
    if zbits < 64:
        ms = [m for m in ms if m % (64 // zbits) == 0]
    return list(ms) or [0]


def nrange(zbits):
    if zbits == 16:
        # x < 2^-16 needs N = 4 xn terms, so only these N are ever
        # dispatched to (xn = 1..5; N = 24 would exceed NMAX)
        return [4, 8, 12, 16, 20]
    return range(1, NMAX + 1)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode.endswith("16"):
        zbits = 16
    elif mode.endswith("32"):
        zbits = 32
    else:
        zbits = 64
    mode = mode.replace("32", "").replace("16", "")
    basename = ("mpn_exp_series" if zbits == 64
                else "mpn_exp_series%d" % zbits)
    allname = ("exp_series_all" if zbits == 64
               else "exp_series%d_all" % zbits)

    out = [header(), emit_tables()]

    if mode == "all":
        entries = []
        for N in nrange(zbits):
            if zbits == 32 and N in HAND32:
                code, eb = HAND32[N]
                out.append(code % "_hand")
                out.append("")
                entries.append(("mpn_exp_series32_%d_hand" % N, N, 1, eb, 0))
            for m in mrange(N, zbits):
                if N <= 3 or m == 0:
                    g = Gen(N, 2 if N >= 4 else 0, zbits,
                            name="%s_%d" % (basename, N))
                    g.m = 0
                    code, eb = g.generate()
                    out.append(code)
                    out.append("")
                    entries.append((g.name, N, 0, eb, 0))
                    continue
                for dr in (0, 1):
                    g = Gen(N, m, zbits, use_divrem=bool(dr))
                    code, eb = g.generate()
                    out.append(code)
                    out.append("")
                    entries.append((g.name, N, m, eb, dr))
        out.append("typedef void (*exp_series_fn)(mp_ptr, mp_srcptr);")
        out.append("typedef struct { exp_series_fn f; int N; int xn; int zbits; int m; int errbound; int dr; }")
        out.append("exp_series_entry;")
        out.append("")
        out.append("exp_series_entry %s[] = {" % allname)
        for (nm, N, m, eb, dr) in entries:
            xn = -(-N * zbits // 64)
            out.append("    { %s, %d, %d, %d, %d, %d, %d }," % (nm, N, xn, zbits, m, eb, dr))
        out.append("};")
        out.append("")
        out.append("int %s_count = %d;" % (allname, len(entries)))
    elif mode == "best":
        with open(sys.argv[2]) as f:
            raw = json.load(f)
        best = {}
        for k, v in raw.items():
            if isinstance(v, dict):
                best[int(k)] = (int(v["m"]), int(v.get("dr", 0)))
            else:
                best[int(k)] = (int(v), 0)
        out.append("typedef void (*exp_series_fn)(mp_ptr, mp_srcptr);")
        out.append("")
        ebs = {}
        n0 = min(nrange(zbits))
        for N in nrange(zbits):
            if zbits == 32 and N in HAND32:
                code, eb = HAND32[N]
                out.append(code % "")
                out.append("")
                ebs[N] = eb
                continue
            m, dr = best.get(N, (0, 0)) if N >= 4 else (0, 0)
            g = Gen(N, m if N >= 4 else 0, zbits, use_divrem=bool(dr),
                    name="%s_%d" % (basename, N))
            if N <= 3:
                g.m = 0
            code, eb = g.generate()
            out.append(code)
            out.append("")
            ebs[N] = eb
        out.append("/* function table indexed by N (%d..%d); xn = ceil(N*%d/64) */"
                   % (n0, NMAX, zbits))
        out.append("const exp_series_fn %s_tab[%d] = {" % (basename, NMAX + 1))
        for N in range(NMAX + 1):
            out.append("    %s," % ("NULL" if N < n0 else "%s_%d" % (basename, N)))
        out.append("};")
        out.append("")
        out.append("/* one-sided ulp error bounds, same indexing */")
        out.append("const int %s_err_tab[%d] = {" % (basename, NMAX + 1))
        for N in range(NMAX + 1):
            out.append("    %d," % (0 if N < n0 else ebs[N]))
        out.append("};")
    else:
        raise SystemExit("mode must be all[32] or best[32]")

    print("\n".join(out))




import re

# ===========================================================================
# Driver for the fixed module: regenerate src/fixed/exp_rs_hard.inc.
#
# The rectangular splitting parameters (m, use_divrem) below were
# selected by benchmarking every candidate on x86-64 (Skylake); to
# retune, generate mode "all", time the variants, and update the
# tables.
# ===========================================================================

WINNERS_64 = {4: (2, 0), 5: (3, 0), 6: (3, 0), 7: (5, 0), 8: (3, 0), 9: (3, 0), 10: (4, 0), 11: (4, 0), 12: (4, 0), 13: (4, 1), 14: (4, 1), 15: (4, 1), 16: (3, 1), 17: (4, 1), 18: (4, 1), 19: (4, 1), 20: (4, 1), 21: (4, 1)}

WINNERS_16 = {4: (4, 0), 8: (4, 0), 12: (4, 0), 16: (4, 0),
              20: (4, 0)}

# Fully specialized series for the small sizes, one per n: the
# reduction parameter r is fixed (and hardcoded into the caller, see
# dev/gen_fixed_exp_bitwise_small.py) and N is the smallest number of
# terms with N r + log2(N!) >= 64 n -- the 1/N! of the last
# coefficient is spent on dropping terms rather than banked as slack.
# n = 1 and n = 2 are NOT generated here: at those sizes the RS
# framework's overhead dwarfs the few word multiplies the series
# needs, so they are hand-written in exp_rs_opt_hand.inc (with N
# reduced to 4 and 7 respectively, plain Horner, sloppy high
# products).  n = 4 is absent because there r = 32 is fastest and its
# N = 2 n is already factorial-minimal.  Entries are n: (r, N, m,
# use_divrem), benchmarked end to end.
WINNERS_OPT = {1: (12, 5, 5, 0), 2: (16, 8, 4, 0), 3: (16, 11, 4, 0),
               4: (16, 14, 4, 0), 5: (16, 17, 4, 0)}

WINNERS_32 = {1: (0, 0), 2: (0, 0), 3: (1, 0), 4: (1, 0), 5: (1, 0), 6: (4, 0), 7: (4, 0), 8: (6, 0), 9: (6, 0), 10: (4, 0), 11: (4, 0), 12: (4, 0), 13: (4, 0), 14: (4, 0), 15: (4, 0), 16: (4, 0), 17: (4, 0), 18: (4, 0), 19: (4, 0), 20: (4, 0), 21: (4, 1)}


def generate_best(zbits, best):
    basename = ("mpn_exp_series" if zbits == 64
                else "mpn_exp_series%d" % zbits)
    out = [header(), emit_tables()]
    out.append("typedef void (*exp_series_fn)(mp_ptr, mp_srcptr);")
    out.append("")
    ebs = {}
    n0 = min(nrange(zbits))
    skipped = set()
    for N in nrange(zbits):
        if zbits == 32 and N % 2 == 1 and N < max(nrange(zbits)):
            # the dispatch reaches the 32-bit-range family only with
            # N = min(2 n, Nmax): odd N below the top entry is dead
            skipped.add(N)
            continue
        if zbits == 32 and N in HAND32:
            code, eb = HAND32[N]
            out.append(code % "")
            out.append("")
            ebs[N] = eb
            continue
        m, dr = best.get(N, (0, 0)) if N >= 4 else (0, 0)
        g = Gen(N, m if N >= 4 else 0, zbits, use_divrem=bool(dr),
                name="%s_%d" % (basename, N))
        if N <= 3:
            g.m = 0
        code, eb = g.generate()
        out.append(code)
        out.append("")
        ebs[N] = eb
    out.append("/* function table indexed by N (%d..%d); xn = ceil(N*%d/64) */"
               % (n0, NMAX, zbits))
    out.append("const exp_series_fn %s_tab[%d] = {" % (basename, NMAX + 1))
    for N in range(NMAX + 1):
        out.append("    %s," % ("NULL" if N not in ebs
                                else "%s_%d" % (basename, N)))
    out.append("};")
    out.append("")
    out.append("/* one-sided ulp error bounds, same indexing */")
    out.append("const int %s_err_tab[%d] = {" % (basename, NMAX + 1))
    for N in range(NMAX + 1):
        out.append("    %d," % (0 if N not in ebs else ebs[N]))
    out.append("};")
    return "\n".join(out)


def _convert_types(s):
    s = s.replace('mp_srcptr', 'nn_srcptr').replace('mp_ptr', 'nn_ptr')
    s = s.replace('mp_limb_t', 'ulong')
    return s


def _drop_err_tabs(s):
    s = re.sub(r'/\* one-sided ulp error bounds[^*]*\*/\n', '', s)
    s = re.sub(r'(static )?const int \w+_err_tab\[\d+\] = \{[^}]*\};\n*', '', s)
    return s


def _statify(s):
    s = re.sub(r'\nvoid\nmpn_exp_series', '\nstatic void\nmpn_exp_series', s)
    s = s.replace('\nconst exp_series_fn ', '\nstatic const exp_series_fn ')
    s = s.replace('typedef void (*exp_series_fn)(mp_ptr, mp_srcptr);', '')
    return s


if __name__ == "__main__":
    # The 2^-16 window family and the per-n opt series died with the
    # restructure: the r >= 32 contract of the public bitwise functions
    # removed the former, and the fully specialized per-size paths in
    # exp_opt_<n>.c carry their own static series (emitted by
    # dev/tune_fixed.py) in place of the latter.  Only the 2^-64 and
    # 2^-32 families remain here.
    f64 = generate_best(64, WINNERS_64)
    f32 = generate_best(32, WINNERS_32)

    f64 = f64[f64.index('/* c_k = 20!/k!'):]
    m32 = re.search(r'mpn_exp_series32_\d+\(', f32)
    i32 = f32.rindex('\n', 0, f32.rindex('void', 0, m32.start()))
    f32 = f32[i32:]

    out = ("/* Generated by dev/gen_fixed_exp_rs_hard.py -- private static\n"
           "   Taylor series routines for the fixed module -- do not edit"
           " by hand. */\n\n"
           "typedef void (*exp_series_fn)(nn_ptr, nn_srcptr);\n\n"
           + _statify(f64) + "\n" + _statify(f32))
    out = _drop_err_tabs(out)
    out = out.replace('mpn_exp_series32_', '_fixed_exp_rs32_')
    out = out.replace('mpn_exp_series_', '_fixed_exp_rs_')
    out = _convert_types(out)
    print(out.rstrip("\n"))
