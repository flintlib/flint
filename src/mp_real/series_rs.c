/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* Tapered rectangular splitting for the reduced series of exp, sin,
   1 - cos, sinh, cosh - 1, atan and atanh at a fraction x < 2^-r of n
   limbs.  With v = x (exp) or v = x^2 and s = r resp. 2r the bits each
   power of v gains,

       P = sum_(k < N) a_k v^k,   a_k = 1/k!, 1/(2k+1)!, 1/(2k+2)!, 1/(2k+1)

   (exp, alternating for exp(-x); sin t / t and sinh t / t; (1 - cos t)/t^2 and (cosh t - 1)/t^2;
   atan t / t and atanh t / t, alternating for the circular functions),
   and the result is P, x P or v P.

   Rectangular splitting in blocks of m (even) terms,

       P = S_0,   S_b = sum_(i < m) a_(mb+i) v^i + v^m S_(b+1),

   with the m powers v, .., v^m computed once and every term a scalar
   multiple c v^i, c a single limb: descending, the factorial
   coefficients chain integrally (a_(k-1)/a_k = F(k) a product of one
   or two small integers), so the sum is carried as U/d with d the
   product of the F since the last division, and one mpn_divrem_1 by d
   closes each run of F that fills a limb; for 1/(2k+1) the terms of a
   group of consecutive odd numbers with product D share the integer
   numerators D/(2k+1), and moving to the next group rescales U by
   D'/D (one mpn_mul_1 and one mpn_divrem_1).

   TAPERING.  Block b enters P multiplied by v^(mb) < 2^(-smb), so it
   runs in a window of the top

       L_b = min(n, ceil((64 n + g - s m b) / 64))

   fraction limbs (plus the units limb), g = bits(N) + 3 guard bits;
   within the block, v^i has floor(s i / 64) leading zero limbs, which
   its addmul_1 skips; the multiplication of S_(b+1) by v^m reads only
   the significant limbs of both and writes the (wider) window of
   block b.  The cost: m full products for the powers, N/m boundary
   products at widths falling linearly to zero (about N/(3m) full
   products), and scalar passes over the shrinking windows -- against
   N/3 full products for the tapered Horner scheme of
   series_tapered.c, which wins only below about a dozen limbs.

   Signs.  In the alternating families U is kept in two's complement
   across the units limb; with m even, every block starts with a
   positive term, so S_b >= 0 at each boundary multiplication, and at a
   division following a negative term, where -d <= U < 0 (the terms
   decrease), U + d is divided and 1 subtracted after (for a rescale,
   U + D and D' subtracted after), keeping every division unsigned.

   Error, in ulps 2^(-64n) of P (then of the result).  An error of delta
   window units committed while the running divisor is d at term k
   reaches P as at most delta a_k / d <= delta window units, multiplied
   by the pending v^(mb) <= 2^(-smb): for b >= 1 each event (at most 3
   units: a window-truncated term, a power's error times a_k <= 1, a
   floor, a boundary product) is below 2^(-min(g, sm)) ulps, and there
   are at most 5 N of them, below 5/8 ulp in all.  In block 0: the
   powers err by at most 1 (v = x^2 by sqrhigh; 0 for v = x) resp. 3
   ulps (the truncated products), weighted by a_i: at most
   3 (e - 2) < 2.2 for exp, 1/6 + 3 (1/5! + ..) < 0.2 for sin and sinh,
   below 0.1 for the cosines and below 1/3 + 3 H < 4.6 for atan and
   atanh with H = 1/5 + 1/7 + .. + 1/(2m - 1) < 1.4 (m <= 32); the
   final division floors at most 1 ulp; an earlier floor reaches P
   divided by the divisor in force at the end of its run, which for
   the factorials means where a_k < 2^-40 (the run filled a limb), and
   for 1/(2k+1) by the next group's product D', above 2^50 except for
   the last group (the groups are formed greedily from the top), so
   one more ulp at most; the last boundary product at most 2 a_m; the
   omitted tail at most 1/8.  So |P - P~| < 4.6 for exp and below
   7.5 for the others.  The final product x P (resp. v P, v with its
   own ulp) damps that by x < 2^-r <= 2^-8 and adds at most 1 ulp for
   the mulhigh (and 1/2 for v's error times P <= 1/2 in the cosines),
   so every result is within 5 ulps for exp and 2 ulps otherwise.
   Outputs may alias x. */

/* a lower bound for 4 log2(F), F >= 1 */
static inline slong
_lg4_lower(ulong F)
{
    static const unsigned char tab[8] = { 0, 0, 1, 1, 2, 2, 2, 3 };
    slong bc = FLINT_BIT_COUNT(F);
    ulong top = (bc >= 4) ? (F >> (bc - 4)) & 7 : (F << (4 - bc)) & 7;
    return 4 * (bc - 1) + tab[top];
}

static inline slong
_rs_window(slong n, slong g, slong sm, slong b)
{
    slong bits = FLINT_BITS * n + g - sm * b;
    slong L = (bits <= 0) ? 1 : (bits + FLINT_BITS - 1) / FLINT_BITS;
    return FLINT_MIN(L, n);
}

/* coefficient families */
#define RS_EXP 0        /* 1/k!                    */
#define RS_ODD 1        /* 1/(2k+1)!               */
#define RS_EVEN 2       /* 1/(2k+2)!               */
#define RS_RECIP 3      /* 1/(2k+1)                */

/* the number of terms: a_N 2^(-sN) < 2^(-64n-3) */
static slong
_rs_terms(int fam, slong n, slong s)
{
    slong N, lg4 = 0, T4 = 4 * (FLINT_BITS * n + 3);

    for (N = 1; ; N++)
    {
        if (fam == RS_EXP)
            lg4 += _lg4_lower((ulong) N);
#if FLINT_BITS == 64
        else if (fam == RS_ODD)
            lg4 += _lg4_lower((ulong) (2 * N + 1) * (ulong) (2 * N));
        else if (fam == RS_EVEN)
            lg4 += _lg4_lower((ulong) (2 * N + 2) * (ulong) (2 * N + 1));
#else
        /* the factors separately (lower bounds add), the product
           possibly beyond a word */
        else if (fam == RS_ODD)
            lg4 += _lg4_lower((ulong) (2 * N + 1)) + _lg4_lower((ulong) (2 * N));
        else if (fam == RS_EVEN)
            lg4 += _lg4_lower((ulong) (2 * N + 2)) + _lg4_lower((ulong) (2 * N + 1));
#endif
        else
            lg4 = _lg4_lower((ulong) (2 * N + 1));
        if (lg4 + 4 * s * N > T4)
            return N;
    }
}

static slong
_rs_blocksize(slong N)
{
    slong m = 2;

    /* m ~ sqrt(N / 1.5), measured within a few percent of the per-size
       optimum from 12 to 256 limbs */
    while (3 * (m + 2) * (m + 2) <= 2 * N + 3 * (m + 1))
        m += 2;
    return FLINT_MIN(m, 32);
}

#define PAD 3

/* P (n + 1 limbs, units limb on top) for the family fam from the powers
   pw(1..m) (pw(i) + j the grid limb j of v^i, PAD zero limbs below) */
static void
_rs_sum(nn_ptr s, nn_srcptr pw, slong pslot, slong n, slong s_bits,
    slong m, int fam, int alt, nn_ptr tmp)
{
    slong N, g, k, b, lo, i, zm, sm;
    ulong d = 1, D = 1, c;
    slong grp_lo = 0;
    int fact = (fam != RS_RECIP);

#define SL(j) (pw + ((j) - 1) * pslot + PAD)

    N = _rs_terms(fam, n, s_bits);
    g = FLINT_BIT_COUNT(N) + 3;
    sm = s_bits * m;
    zm = sm / FLINT_BITS;

    flint_mpn_zero(s, n + 2);
    b = (N - 1) / m;
    lo = n - _rs_window(n, g, sm, b);

    for (k = N - 1; k >= 0; k--)
    {
        int neg = alt && (k & 1);
        int pneg = alt && !(k & 1);     /* the previous term, k + 1 */

        i = k - m * b;

        if (fact)
        {
            if (k < N - 1)
            {
                ulong F, hi, lw;

                /* F fits a word: 2k + 4 < 2^(FLINT_BITS/2), which the
                   term counts keep far from (k ~ 64n/s, s >= 8) */
                FLINT_ASSERT(FLINT_BITS == 64 || 2 * k + 4 < (WORD(1) << 16));
                if (fam == RS_EXP)
                    F = (ulong) (k + 1);
                else if (fam == RS_ODD)
                    F = (ulong) (2 * k + 3) * (ulong) (2 * k + 2);
                else
                    F = (ulong) (2 * k + 4) * (ulong) (2 * k + 3);

                umul_ppmm(hi, lw, d, F);
                if (hi != 0)
                {
                    if (pneg)
                        s[n] += d;
                    mpn_divrem_1(s + lo, 0, s + lo, n - lo + 1, d);
                    if (pneg)
                        s[n] -= 1;
                    d = F;
                }
                else
                    d = lw;
            }
            c = d;
        }
        else
        {
            if (k == N - 1 || k < grp_lo)
            {
                ulong ND = 1, hi, lw;
                slong j = k;

                while (j >= 0)
                {
                    umul_ppmm(hi, lw, ND, (ulong) (2 * j + 1));
                    if (hi != 0)
                        break;
                    ND = lw;
                    j--;
                }
                grp_lo = j + 1;

                if (k < N - 1)
                {
                    if (pneg)
                        s[n] += D;
                    s[n + 1] = mpn_mul_1(s + lo, s + lo, n - lo + 1, ND);
                    mpn_divrem_1(s + lo, 0, s + lo, n - lo + 2, D);
                    if (pneg)
                        s[n] -= ND;
                }
                D = ND;
            }
            c = D / (ulong) (2 * k + 1);
        }

        if (i == 0)
        {
            if (neg)
                s[n] -= c;
            else
                s[n] += c;
        }
        else
        {
            slong zi = (s_bits * i) / FLINT_BITS;
            slong len = n - zi - lo;

            if (len > 0)
            {
                ulong cy;
                slong p = n - zi;

                if (neg)
                {
                    cy = mpn_submul_1(s + lo, SL(i) + lo, len, c);
                    while (cy != 0 && p <= n)
                    {
                        ulong t = s[p];
                        s[p] = t - cy;
                        cy = (t < cy);
                        p++;
                    }
                }
                else
                {
                    cy = mpn_addmul_1(s + lo, SL(i) + lo, len, c);
                    while (cy != 0 && p <= n)
                    {
                        ulong t = s[p] + cy;
                        cy = (t < cy);
                        s[p] = t;
                        p++;
                    }
                }
            }
        }

        if (i == 0 && k > 0)
        {
            /* S_b = s v^m on the window of block b - 1: the product of
               the top W limbs of s (zero below its window) and of v^m's
               significant part lands at grid position n - (W - 1) - zm,
               at or below the new window's bottom */
            slong Lo = n - lo, Ln, W, base;

            b--;
            Ln = _rs_window(n, g, sm, b);
            lo = n - Ln;
            if (zm > Ln)
            {
                /* v^m vanishes on the window */
                flint_mpn_zero(s + lo, Ln + 1);
            }
            else
            {
                /* Lo + zm <= n + 2 as b + 1 >= 1, so the read of v^m
                   stays within its PAD zero limbs */
                W = FLINT_MAX(Lo, Ln - zm) + 1;
                flint_mpn_mulhigh_n(tmp, s + n - W + 1, SL(m) + n - zm - W, W);
                base = n - (W - 1) - zm;
                flint_mpn_copyi(s + lo, tmp + (lo - base), n - zm - lo + 1);
                flint_mpn_zero(s + n - zm + 1, zm);
            }
        }
    }

    mpn_divrem_1(s, 0, s, n + 1, fact ? d : D);
    if (fam == RS_EVEN)
        mpn_rshift(s, s, n + 1, 1);     /* a_0 = 1/2 */

#undef SL
}

/* the powers v^2 .. v^m of v = pw(1) (grid limbs, pw(1) already set) */
static void
_rs_powers(nn_ptr pw, slong pslot, slong n, slong s_bits, slong m, slong pad)
{
    slong i;

#define SL(j) (pw + ((j) - 1) * pslot + pad)
    for (i = 1; i <= m; i++)
        flint_mpn_zero(SL(i) - pad, pad);

    for (i = 2; i <= m; i++)
    {
        slong a = i / 2, bb = i - a;
        slong za = (s_bits * a) / FLINT_BITS, zb = (s_bits * bb) / FLINT_BITS;
        slong K = n - za - zb;

        /* v^a v^b from the significant limbs of each: floor(s a / 64)
           resp. floor(s b / 64) top limbs are zero */
        if (K <= 0)
            flint_mpn_zero(SL(i), n);
        else
        {
            if (a == bb)
                flint_mpn_sqrhigh(SL(i), SL(a) + zb, K);
            else
                flint_mpn_mulhigh_n(SL(i), SL(a) + zb, SL(bb) + za, K);
            flint_mpn_zero(SL(i) + K, n - K);
        }
    }
#undef SL
}

void
_mp_real_series_rs_sin_cos(nn_ptr ysin, nn_ptr yg, nn_srcptr x, slong n,
    flint_bitcnt_t r, int hyperbolic)
{
    slong s_bits = 2 * (slong) r, N, m, pslot;
    nn_ptr pw, s, tmp, z;
    int alt = !hyperbolic;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r >= 8);

    N = _rs_terms(RS_ODD, n, s_bits);
    m = _rs_blocksize(N);

    TMP_START;
    pslot = n + PAD;
    pw = TMP_ALLOC((m * pslot + 3 * (n + 3)) * sizeof(ulong));
    s = pw + m * pslot;
    tmp = s + 2 * (n + 3);
    z = pw + PAD;

    flint_mpn_sqrhigh(z, x, n);
    _rs_powers(pw, pslot, n, s_bits, m, PAD);

    /* the products go through tmp: ysin and yg may alias x */
    if (yg != NULL)
    {
        _rs_sum(s, pw, pslot, n, s_bits, m, RS_EVEN, alt, tmp);
        flint_mpn_mulhigh_n(tmp, z, s, n);  /* s < 1 */
        if (ysin == NULL)
        {
            flint_mpn_copyi(yg, tmp, n);
            TMP_END;
            return;
        }
        flint_mpn_copyi(s + n + 3, tmp, n);
    }

    if (ysin != NULL)
    {
        _rs_sum(s, pw, pslot, n, s_bits, m, RS_ODD, alt, tmp);
        flint_mpn_mulhigh_n(tmp, x, s, n);
        if (s[n] != 0)
            mpn_addmul_1(tmp, x, n, s[n]);
        flint_mpn_copyi(ysin, tmp, n);
        if (yg != NULL)
            flint_mpn_copyi(yg, s + n + 3, n);
    }

    TMP_END;
}

void
_mp_real_series_rs(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r,
    int func)
{
    slong s_bits, N, m, pslot;
    int fam, alt;
    nn_ptr pw, s, tmp;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r >= 8);

    if (func == MP_REAL_SERIES_SIN || func == MP_REAL_SERIES_SINH)
    {
        _mp_real_series_rs_sin_cos(res, NULL, x, n, r,
            func == MP_REAL_SERIES_SINH);
        return;
    }

    if (func == MP_REAL_SERIES_COS || func == MP_REAL_SERIES_COSH)
    {
        _mp_real_series_rs_sin_cos(NULL, res, x, n, r,
            func == MP_REAL_SERIES_COSH);
        return;
    }

    if (func == MP_REAL_SERIES_EXP || func == MP_REAL_SERIES_EXP_NEG)
    {
        fam = RS_EXP;
        alt = (func == MP_REAL_SERIES_EXP_NEG);
        s_bits = r;
    }
    else
    {
        FLINT_ASSERT(func == MP_REAL_SERIES_ATAN || func == MP_REAL_SERIES_ATANH);
        fam = RS_RECIP;
        alt = (func == MP_REAL_SERIES_ATAN);
        s_bits = 2 * (slong) r;
    }

    N = _rs_terms(fam, n, s_bits);

    /* exp x = 1 + x resp. atan(h) x = x to the working precision */
    if (N <= (fam == RS_EXP ? 2 : 1))
    {
        if (fam == RS_EXP && alt)
        {
            /* exp(-x) = 1 - x */
            res[n] = flint_mpn_zero_p(x, n);
            mpn_neg(res, x, n);
            return;
        }
        flint_mpn_copyi(res, x, n);
        if (fam == RS_EXP)
            res[n] = 1;
        return;
    }

    m = _rs_blocksize(N);

    TMP_START;
    pslot = n + PAD;
    pw = TMP_ALLOC((m * pslot + 2 * (n + 3)) * sizeof(ulong));
    s = pw + m * pslot;
    tmp = s + (n + 3);

    if (fam == RS_EXP)
        flint_mpn_copyi(pw + PAD, x, n);
    else
        flint_mpn_sqrhigh(pw + PAD, x, n);
    _rs_powers(pw, pslot, n, s_bits, m, PAD);

    _rs_sum(s, pw, pslot, n, s_bits, m, fam, alt, tmp);

    if (fam == RS_EXP)
        flint_mpn_copyi(res, s, n + 1);
    else
    {
        /* through tmp: res may alias x */
        flint_mpn_mulhigh_n(tmp, x, s, n);
        if (s[n] != 0)
            mpn_addmul_1(tmp, x, n, s[n]);
        flint_mpn_copyi(res, tmp, n);
    }

    TMP_END;
}

/* ==== the tangent ==========================================================

   tan t = t + t^3 S, S = sum_k c_k z^k, z = t^2, c_k = T_(k+2)/(2k+3)!.
   The coefficients are not ratios of small integers, but the leading
   ones share small denominators: with Q_j the lcm of the reduced
   denominators of c_0 .. c_k, chunk j (elem_tables.c) holds the terms
   whose Q_j fits j + 1 limbs (12 terms in one limb, then 5-7 terms per
   additional limb on 64-bit machines), with the integer numerators
   N_k = c_k Q_j.  One rectangular splitting runs over the chunked terms
   k < K (blocks of m, K a multiple of m when a tail follows),
   descending, carrying Q_j times the partial sum with j + 1 units
   limbs: each term is j + 1 addmul_1 passes, crossing into chunk j - 1
   one exact-in-the-units division by R_j = Q_j / Q_(j-1), and the end
   one division by Q_0.  A chunk is used while its width is at most a
   third of the window at its first term (where the scalar passes stop
   being cheaper than products); the terms from K on are the tapered
   Horner tail of series_tapered.c (whose tables end at
   MP_REAL_SERIES_TAN_NMAX limbs: beyond, the chunks must cover every
   term, see _mp_real_series_rs_tan_ok), entering as the top block's
   continuation Q_j V_K.  Windows as in _mp_real_series_rs, for the
   term k at t^3 z^k < 2^(-r(2k+3)); all errors in S are damped by
   t^3, so the result is within 3 ulps (the products t^3 = z t and
   t^3 S). */

#define TAN_PAD (2 * MP_REAL_SERIES_TAN_RS_J + 3)

static inline slong
_tan_window(slong n, slong g, slong r, slong k)
{
    slong bits = FLINT_BITS * n + g - r * (2 * k + 3);
    slong L = (bits <= 0) ? 1 : (bits + FLINT_BITS - 1) / FLINT_BITS;
    return FLINT_MIN(L, n);
}

/* (s, grid [lo_new, n + w)) = z^m times (s, grid [lo_old, n + w)) */
static void
_tan_boundary(nn_ptr s, slong n, slong w, slong lo_old, slong Ln,
    nn_srcptr zmp, slong zm, nn_ptr tmp)
{
    slong lo_new = n - Ln, W, base, top;

    if (zm >= Ln + w)
    {
        flint_mpn_zero(s + lo_new, Ln + w);
        return;
    }

    W = FLINT_MAX(n + w - lo_old, Ln + w - zm);
    flint_mpn_mulhigh_n(tmp, s + n + w - W, zmp + n - zm - W, W);
    base = n + w - zm - W;
    top = n + w - zm - 1;
    flint_mpn_copyi(s + lo_new, tmp + (lo_new - base), top - lo_new + 1);
    flint_mpn_zero(s + top + 1, zm);
}

/* the number of terms from the Horner table's bounds for log2 c_k: the
   first omitted term below 2^(-64 n - 2); MP_REAL_SERIES_TAN_K if more */
static slong
_tan_terms(slong n, slong r)
{
    slong N;

    for (N = 0; N < MP_REAL_SERIES_TAN_K; N++)
        if (_mp_real_series_tan_lg[N] - r * (3 + 2 * N) < -FLINT_BITS * n - 2)
            break;
    return FLINT_MAX(N, 1);
}

int
_mp_real_series_rs_tan_ok(slong n, flint_bitcnt_t r)
{
    slong N;

    if ((slong) r < MP_REAL_SERIES_TAN_RMIN)
        return 0;
    if (n <= MP_REAL_SERIES_TAN_NMAX)
        return 1;
    N = _tan_terms(n, (slong) r);
    return N < MP_REAL_SERIES_TAN_K
        && N <= _mp_real_series_tan_rs_start[MP_REAL_SERIES_TAN_RS_J];
}

void
_mp_real_series_rs_tan(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r)
{
    const short * start = _mp_real_series_tan_rs_start;
    slong N, g, J, K, Kc, m, zm, pslot, b, lo, k, i, w, l, jt;
    slong rr = (slong) r;
    nn_ptr pw, s, tmp, z;
    TMP_INIT;

    FLINT_ASSERT(n >= 1 && _mp_real_series_rs_tan_ok(n, r));

    /* the Horner scheme below 10 limbs, and below 20 at r >= 64 (few
       terms): measured 1.0-1.15 times faster there, the chunks 1.05-1.9
       times faster above (the width rule and the block size flat within
       a few percent) */
    if (n < 10 || (n < 20 && rr >= 64))
    {
        _mp_real_series_tapered(res, x, n, r, MP_REAL_SERIES_TAN);
        return;
    }

    N = _tan_terms(n, rr);
    g = FLINT_BIT_COUNT(N) + 3;

    J = 0;
    if (n <= MP_REAL_SERIES_TAN_NMAX)
    {
        while (J < MP_REAL_SERIES_TAN_RS_J && start[J] < N
                && 3 * (J + 1) <= _tan_window(n, g, rr, start[J]))
            J++;
    }
    else
    {
        /* no Horner table at this size: the chunks cover every term */
        while (start[J] < N)
            J++;
    }

    Kc = (J == 0) ? 0 : FLINT_MIN(N, start[J]);
    m = _rs_blocksize(Kc);
    K = (Kc < N) ? (Kc / m) * m : N;

    if (K < 2)
    {
        _mp_real_series_tapered(res, x, n, r, MP_REAL_SERIES_TAN);
        return;
    }

    /* the chunk of the top term */
    for (jt = 0; start[jt + 1] <= K - 1; jt++)
        ;
    w = jt + 1;

    TMP_START;
    pslot = n + TAN_PAD;
    pw = TMP_ALLOC((m * pslot + 3 * (n + TAN_PAD)) * sizeof(ulong));
    s = pw + m * pslot;
    tmp = s + (n + TAN_PAD);
    z = pw + TAN_PAD;
    zm = (2 * rr * m) / FLINT_BITS;

#define SL(j) (pw + ((j) - 1) * pslot + TAN_PAD)

    flint_mpn_sqrhigh(z, x, n);
    _rs_powers(pw, pslot, n, 2 * rr, m, TAN_PAD);
    flint_mpn_zero(s, n + TAN_PAD);

    b = (K - 1) / m;
    if (K < N)
    {
        /* S_(b+1) = V_K, the Horner tail, times Q_jt */
        nn_srcptr q = _mp_real_series_tan_rs_q + jt * (jt + 1) / 2;
        slong cur = _mp_real_series_tapered_horner(tmp, z, n, r,
            MP_REAL_SERIES_TAN, K, N, g);

        lo = n - cur;
        for (l = 0; l < w; l++)
        {
            ulong cy = mpn_addmul_1(s + lo + l, tmp, cur, q[l]);
            slong p = n + l;
            while (cy != 0 && p < n + w)
            {
                ulong t = s[p] + cy;
                cy = (t < cy);
                s[p] = t;
                p++;
            }
        }
        _tan_boundary(s, n, w, lo, _tan_window(n, g, rr, m * b), SL(m), zm, tmp);
    }
    lo = n - _tan_window(n, g, rr, m * b);

    for (k = K - 1; k >= 0; k--)
    {
        nn_srcptr num;

        i = k - m * b;

        if (k < start[w - 1])
        {
            /* into chunk w - 2: divide by R_(w-1) = Q_(w-1) / Q_(w-2) */
            nn_srcptr R = _mp_real_series_tan_rs_r + (w - 1) * w / 2;
            slong rn = w, len = n + w - lo;

            while (rn > 1 && R[rn - 1] == 0)
                rn--;
            if (rn == 1)
                mpn_divrem_1(s + lo, 0, s + lo, len, R[0]);
            else
            {
                mpn_tdiv_qr(tmp, tmp + len, 0, s + lo, len, R, rn);
                flint_mpn_copyi(s + lo, tmp, len - rn + 1);
                flint_mpn_zero(s + lo + len - rn + 1, rn - 1);
            }
            w--;
        }

        num = _mp_real_series_tan_rs_num + _mp_real_series_tan_rs_num_off[k];

        if (i == 0)
            mpn_add_n(s + n, s + n, num, w);
        else
        {
            slong zi = (2 * rr * i) / FLINT_BITS;

            for (l = 0; l < w; l++)
            {
                slong len = n - zi - lo + l, p;
                ulong cy;

                if (num[l] == 0 || len <= 0)
                    continue;
                /* num[l] B^l z^i: z^i's limbs from lo - l land from lo */
                cy = mpn_addmul_1(s + lo, SL(i) + lo - l, len, num[l]);
                p = lo + len;
                while (cy != 0 && p < n + w)
                {
                    ulong t = s[p] + cy;
                    cy = (t < cy);
                    s[p] = t;
                    p++;
                }
            }
        }

        if (i == 0 && k > 0)
        {
            slong lo_old = lo;
            b--;
            lo = n - _tan_window(n, g, rr, m * b);
            _tan_boundary(s, n, w, lo_old, n - lo, SL(m), zm, tmp);
        }
    }

    FLINT_ASSERT(w == 1);

    /* S = U / Q_0 < 1/2, then tan t = t + (z t) S */
    mpn_divrem_1(s, 0, s, n + 1, _mp_real_series_tan_rs_q[0]);
    FLINT_ASSERT(s[n] == 0);
    flint_mpn_mulhigh_n(tmp, z, x, n);
    flint_mpn_mulhigh_n(tmp + n, tmp, s, n);
    mpn_add_n(res, x, tmp + n, n);

#undef SL
    TMP_END;
}
