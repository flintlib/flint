/*
    Copyright (C) 2014, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

/* Pure mpn port of _arb_exp_sum_bs_powtab (single-threaded):

       T / (Q 2^Qexp) = sum_{k=1}^{N} x^k / (k! 2^(r k)),

   x > 0 given as (xp, xn) limbs, by binary splitting with a table of
   the powers x^step that the fixed midpoint recursion needs (the
   b - a <= 2 cases inlined), the 2-adic parts of the factorials
   folded into Qexp.  All numbers are nonnegative mpn integers with
   tracked lengths; every flint_mpn_mul call orders its operands.

   The split-step exponent computation and the power-table structure
   are verbatim from the arb original, and t-exp_sum_bs checks the
   port against it for bit identity. */

/* exponents (m - a) reachable from [0, n) with m = a + (b-a)/2 and
   b - a = 2 inlined; at most 2 FLINT_BITS entries, returned sorted */
static slong
compute_bs_exponents(slong * tab, slong n)
{
    slong a, b, aa, ba, bb, length;

    if (n == 1)
    {
        tab[0] = 1;
        return 1;
    }
    if (n == 2 || n == 3 || n == 4)
    {
        tab[0] = 1;
        tab[1] = 2;
        return 2;
    }
    if (n == 6)
    {
        tab[0] = 1;
        tab[1] = 2;
        tab[2] = 3;
        return 3;
    }

    a = n >> 1;
    b = n - (n >> 1);
    tab[0] = a;
    length = 1;

    for (;;)
    {
        aa = a >> 1;
        ba = b >> 1;
        bb = b - ba;

        tab[length] = ba;
        length++;

        if (ba == 3)
        {
            tab[length] = 2;
            tab[length + 1] = 1;
            length += 2;
            break;
        }

        if (ba == 1 || (ba == 2 && (n & (n - 1)) == 0))
            break;

        if (aa != ba && aa != 1)
        {
            tab[length] = aa;
            length++;
        }

        a = aa;
        b = bb;
    }

    if (tab[length - 1] != 1)
    {
        tab[length] = 1;
        length++;
    }

    for (a = 0; a < length / 2; a++)
    {
        b = tab[a];
        tab[a] = tab[length - a - 1];
        tab[length - a - 1] = b;
    }

    return length;
}

static slong
get_exp_pos(const slong * tab, slong step)
{
    slong i;
    for (i = 0; ; i++)
    {
        if (tab[i] == step)
            return i;
        FLINT_ASSERT(tab[i] != 0);
    }
}


/* Number of series terms for the binary splitting at |x| < 2^-r and
   target precision prec bits, accounting for the factorial growth
   (not only the geometric decay -- decisive for small r at high
   precision), then padded to a higher 2-valuation so that the fixed
   midpoint splitting needs fewer distinct entries in the power
   table.  Port of bs_num_terms / _arb_exp_taylor_bound; the short
   table bounds lg(1/N!) from above, and beyond it a floor-log sum
   continues the bound safely (each lg j >= floor(lg j)). */

#define FAC_TABSIZE 256

static const short rec_fac_bound_tab[FAC_TABSIZE] =
{
    0, 0, -1, -2, -4, -6, -9, -12, -15, -18, -21, -25, -28, -32, -36,
    -40, -44, -48, -52, -56, -61, -65, -69, -74, -79, -83, -88, -93, -97,
    -102, -107, -112, -117, -122, -127, -132, -138, -143, -148, -153,
    -159, -164, -169, -175, -180, -186, -191, -197, -202, -208, -214,
    -219, -225, -231, -237, -242, -248, -254, -260, -266, -272, -278,
    -284, -289, -295, -302, -308, -314, -320, -326, -332, -338, -344,
    -350, -357, -363, -369, -375, -382, -388, -394, -401, -407, -413,
    -420, -426, -433, -439, -446, -452, -458, -465, -472, -478, -485,
    -491, -498, -504, -511, -518, -524, -531, -538, -544, -551, -558,
    -564, -571, -578, -585, -591, -598, -605, -612, -619, -626, -632,
    -639, -646, -653, -660, -667, -674, -681, -688, -695, -702, -709,
    -716, -723, -730, -737, -744, -751, -758, -765, -772, -779, -786,
    -793, -801, -808, -815, -822, -829, -836, -844, -851, -858, -865,
    -872, -880, -887, -894, -901, -909, -916, -923, -931, -938, -945,
    -952, -960, -967, -975, -982, -989, -997, -1004, -1011, -1019, -1026,
    -1034, -1041, -1049, -1056, -1064, -1071, -1078, -1086, -1093, -1101,
    -1108, -1116, -1123, -1131, -1139, -1146, -1154, -1161, -1169, -1176,
    -1184, -1192, -1199, -1207, -1214, -1222, -1230, -1237, -1245, -1253,
    -1260, -1268, -1276, -1283, -1291, -1299, -1306, -1314, -1322, -1329,
    -1337, -1345, -1353, -1360, -1368, -1376, -1384, -1391, -1399, -1407,
    -1415, -1423, -1430, -1438, -1446, -1454, -1462, -1470, -1477, -1485,
    -1493, -1501, -1509, -1517, -1525, -1532, -1540, -1548, -1556, -1564,
    -1572, -1580, -1588, -1596, -1604, -1612, -1620, -1628, -1636, -1644,
    -1652, -1660, -1668, -1675
};

slong
_fixed_exp_bs_num_terms(flint_bitcnt_t r, slong prec)
{
    slong N, fb;

    fb = 0;
    for (N = 1; ; N++)
    {
        if (N < FAC_TABSIZE)
            fb = rec_fac_bound_tab[N];
        else
            fb -= FLINT_BIT_COUNT((ulong) N) - 1;

        if (-(slong) r * N + fb < -prec - 1)
            break;
    }

    if (N > 10000)
        while (N % 128 != 0)
            N++;
    if (N > 1000)
        while (N % 16 != 0)
            N++;
    if (N > 100)
        while (N % 2 != 0)
            N++;

    return N;
}

/* ---- small mpn-number helpers: (ptr, len), len 0 means zero ---- */

static slong
nnn_normalize(nn_srcptr p, slong n)
{
    while (n > 0 && p[n - 1] == 0)
        n--;
    return n;
}

/* z = a * b, ordered; z not aliased; returns length */
static slong
nnn_mul(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn)
{
    if (an == 0 || bn == 0)
        return 0;
    if (an >= bn)
        flint_mpn_mul(z, a, an, b, bn);
    else
        flint_mpn_mul(z, b, bn, a, an);
    return nnn_normalize(z, an + bn);
}

/* strip trailing zero LIMBS of q into the limb exponent *e: the
   limb-granular 2-valuation normalization of the base-B kernel
   (the factorial product Q = prod k carries val_2(N!) ~ N bits of
   2-adic content, i.e. whole dead limbs once N exceeds ~64, which
   would otherwise ride along through every Q Q2 and T Q2 product) */
static slong
nnn_strip_low(nn_ptr q, slong qn, slong * e)
{
    slong z = 0;
    while (z < qn - 1 && q[z] == 0)
        z++;
    if (z > 0)
    {
        flint_mpn_copyi(q, q + z, qn - z);
        *e += z;
        qn -= z;
    }
    return qn;
}

/* power table: xpow[i] = x^(xexp[i]), by squarings, preferring
   squares of the previous entries (same structure as the arb
   original, which guarantees one of the four cases applies).
   Returns the malloc'd backing storage (caller frees).

   When xoff is non-NULL, each entry is stored with its trailing
   zero LIMBS stripped and xpow[i] B^xoff[i] = x^xexp[i]: the
   2-adic content of x multiplies under powering
   (val_2(x^e) = e val_2(x)), so even a sub-limb val_2(x) turns
   into whole dead limbs at higher powers, and squaring a stripped
   entry inherits twice its offset before fresh stripping.  With
   xoff NULL nothing is stripped (the bit-granular kernel keeps
   its arb-identical layout). */
static nn_ptr
build_xpow_table(nn_srcptr xp, slong xn, const slong * xexp,
    slong length, nn_srcptr * xpow, slong * xlen, slong * xoff)
{
    nn_ptr storage;
    slong i, off;

    {
        slong total = 0;
        for (i = 0; i < length; i++)
            total += xexp[i] * xn + 1;
        storage = flint_malloc(total * sizeof(ulong));
    }
    off = 0;
    for (i = 0; i < length; i++)
    {
        nn_ptr dst = storage + off;
        off += xexp[i] * xn + 1;

        slong off2 = 0;

        if (i == 0)
        {
            FLINT_ASSERT(xexp[0] == 1);
            flint_mpn_copyi(dst, xp, xn);
            xlen[0] = xn;
        }
        else if (xexp[i] == 2 * xexp[i - 1])
        {
            flint_mpn_sqr(dst, xpow[i - 1], xlen[i - 1]);
            xlen[i] = nnn_normalize(dst, 2 * xlen[i - 1]);
            if (xoff)
                off2 = 2 * xoff[i - 1];
        }
        else if (i >= 2 && xexp[i] == 2 * xexp[i - 2])
        {
            flint_mpn_sqr(dst, xpow[i - 2], xlen[i - 2]);
            xlen[i] = nnn_normalize(dst, 2 * xlen[i - 2]);
            if (xoff)
                off2 = 2 * xoff[i - 2];
        }
        else
        {
            /* an odd 2k+1 case: square then multiply by x, via a
               scratch for the square */
            nn_ptr sq;
            slong sl, src;
            if (xexp[i] == 2 * xexp[i - 1] + 1)
                src = i - 1;
            else
            {
                FLINT_ASSERT(i >= 2 && xexp[i] == 2 * xexp[i - 2] + 1);
                src = i - 2;
            }
            sq = flint_malloc((2 * xlen[src] + 1) * sizeof(ulong));
            flint_mpn_sqr(sq, xpow[src], xlen[src]);
            sl = nnn_normalize(sq, 2 * xlen[src]);
            xlen[i] = nnn_mul(dst, sq, sl, xp, xn);
            flint_free(sq);
            if (xoff)
                off2 = 2 * xoff[src] + xoff[0];
        }
        if (xoff)
        {
            xoff[i] = off2;
            xlen[i] = nnn_strip_low(dst, xlen[i], &xoff[i]);
        }
        xpow[i] = dst;
    }
    return storage;
}

/* ---- the recursion ---- */

typedef struct
{
    const slong * xexp;
    nn_srcptr * xpow;
    const slong * xlen;
    const slong * xoff;     /* per-entry limb offsets when the
                               power table is stored stripped
                               (base-B kernel); NULL otherwise */
    flint_bitcnt_t r;
    slong xn;               /* limbs of x, for size bounds */
}
bs_args;

static slong
qbound(const bs_args * args, slong a, slong b)
{
    return ((b - a) * FLINT_BIT_COUNT((ulong) b + 1)) / FLINT_BITS + 3;
}

/* ==== base-2^FLINT_BITS twin ============================================

   Same sum, but the common power-of-two exponent QE is tracked in
   LIMBS rather than bits and the argument frame is a whole number of
   limbs:

       T / (Q B^QE) = sum_{k=1}^{N} x^k / (k! B^(D k)),   B = 2^64,

   x = (xp, xn) with value x B^-D, D >= 1 limbs.  All 2-valuation
   normalization happens at limb granularity: the factorial factors
   keep their 2-adic content inside Q (no ctz stripping -- Q is at
   most the same worst-case size, prod k <= N^(range), merely up to
   ~N bits denser in practice), and every << Qexp in the recursion
   becomes a zero-and-copy at a limb offset, with no mpn_lshift
   anywhere.  The trade, known from the dev notes on limb-aligned
   burst slicing: the splitting integers scale with the FRAME 64 D
   instead of the argument's true rate r, a content factor
   (64 D) / r that matters for r well below the limb size and fades
   as r grows. */

/* z (length zn) += a B^q, i.e. the addition itself is performed at
   the limb offset q -- no zeroing or copying of shifted operands
   except to fill a gap when z does not reach the offset; returns
   the new length */
static slong
nnn_add_shifted(nn_ptr z, slong zn, nn_srcptr a, slong an, slong q)
{
    ulong cy;

    if (an == 0)
        return zn;
    if (zn <= q)
    {
        flint_mpn_zero(z + zn, q - zn);
        flint_mpn_copyi(z + q, a, an);
        return q + an;
    }
    if (zn - q >= an)
    {
        cy = mpn_add(z + q, z + q, zn - q, a, an);
        z[zn] = cy;
        return zn + (cy != 0);
    }
    cy = mpn_add(z + q, a, an, z + q, zn - q);
    z[q + an] = cy;
    return q + an + (cy != 0);
}

/* T bound in limbs for the range [a, b): T < Q B^QE with
   QE = (b - a) D and Q < b^(b-a) */
static slong
tbound(const bs_args * args, slong D, slong a, slong b)
{
    (void) args;
    return (b - a) * (D + 2) + 4;
}

static void
bsplit(nn_ptr T, slong * tn, nn_ptr Q, slong * qn, slong * QE,
    const bs_args * args, slong D, slong a, slong b)
{
    if (b - a == 1)
    {
        Q[0] = (ulong) (a + 1);
        *qn = 1;
        *QE = D;
        flint_mpn_zero(T, args->xoff[0]);
        flint_mpn_copyi(T + args->xoff[0], args->xpow[0],
            args->xlen[0]);
        *tn = args->xoff[0] + args->xlen[0];
    }
    else if (b - a == 2)
    {
        /* T = x (a+2) B^D + x^2: the low term x^2 goes into T
           directly and the high term is ADDED at limb offset D */
        nn_ptr t1;
        slong l1;
        ulong hi;
        TMP_INIT;
        TMP_START;
        t1 = TMP_ALLOC((args->xlen[0] + 2) * sizeof(ulong));
        t1[args->xlen[0]] = mpn_mul_1(t1, args->xpow[0],
            args->xlen[0], (ulong) (a + 2));
        l1 = nnn_normalize(t1, args->xlen[0] + 1);
        FLINT_ASSERT(args->xexp[1] == 2);
        flint_mpn_zero(T, args->xoff[1]);
        flint_mpn_copyi(T + args->xoff[1], args->xpow[1],
            args->xlen[1]);
        *tn = nnn_add_shifted(T, args->xoff[1] + args->xlen[1],
            t1, l1, D + args->xoff[0]);
        TMP_END;

        umul_ppmm(hi, Q[0], (ulong) (a + 1), (ulong) (a + 2));
        Q[1] = hi;
        *qn = 1 + (hi != 0);
        *QE = 2 * D;
    }
    else
    {
        slong step, m, i, t2n, q2n, l, Q2E;
        nn_ptr T2, Q2, sc;
        TMP_INIT;

        step = (b - a) / 2;
        m = a + step;

        TMP_START;
        T2 = TMP_ALLOC((tbound(args, D, m, b)
            + qbound(args, m, b)
            + tbound(args, D, a, b) + 8) * sizeof(ulong));
        Q2 = T2 + tbound(args, D, m, b);
        sc = Q2 + qbound(args, m, b);

        bsplit(T, tn, Q, qn, QE, args, D, a, m);
        bsplit(T2, &t2n, Q2, &q2n, &Q2E, args, D, m, b);

        /* T = (T Q2) B^Q2E + x^step T2: the high product goes to
           scratch (it reads T), the low product x^step T2 then
           lands in T directly, and the high part is ADDED at the
           limb offset Q2E -- no shifted zero-and-copy */
        l = nnn_mul(sc, T, *tn, Q2, q2n);
        i = get_exp_pos(args->xexp, step);
        flint_mpn_zero(T, args->xoff[i]);
        *tn = args->xoff[i] + nnn_mul(T + args->xoff[i],
            args->xpow[i], args->xlen[i], T2, t2n);
        *tn = nnn_add_shifted(T, *tn, sc, l, Q2E);

        /* Q = Q Q2, trailing zero limbs stripped into QE */
        l = nnn_mul(sc, Q, *qn, Q2, q2n);
        flint_mpn_copyi(Q, sc, l);
        *QE = *QE + Q2E;
        *qn = nnn_strip_low(Q, l, QE);
        TMP_END;
    }
}

/* T (*tn) / (Q (*qn) B^QE) = sum_{k=1}^N x^k / (k! B^(D k)).
   T must have room for N (D + 2) + 4 limbs and Q for
   ((N bits(N+1)) / 64 + 3) limbs. */
void
_fixed_exp_sum_bs_powtab(nn_ptr T, slong * tn, nn_ptr Q,
    slong * qn, slong * QE, nn_srcptr xp, slong xn, slong D,
    slong N)
{
    slong xexp[2 * FLINT_BITS];
    slong xlen[2 * FLINT_BITS];
    slong xoff[2 * FLINT_BITS];
    nn_srcptr xpow[2 * FLINT_BITS];
    nn_ptr storage;
    slong length, i;
    bs_args args;

    FLINT_ASSERT(N >= 1 && xn >= 1 && xp[xn - 1] != 0);
    FLINT_ASSERT(xn <= D);      /* value x B^-D below 1 */

    for (i = 0; i < 2 * FLINT_BITS; i++)
        xexp[i] = 0;
    length = compute_bs_exponents(xexp, N);

    storage = build_xpow_table(xp, xn, xexp, length, xpow, xlen,
        xoff);

    args.xexp = xexp;
    args.xpow = xpow;
    args.xlen = xlen;
    args.xoff = xoff;
    args.r = 0;
    args.xn = xn;

    bsplit(T, tn, Q, qn, QE, &args, D, 0, N);

    flint_free(storage);
}
