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
#include "fixed.h"

/* Rectangular splitting evaluation of sin, cos, sinh, cosh, atan and
   atanh on [0, 2^-32) resp. [0, 2^-64).  All eight entry points share
   the machinery below (the series substitute z = x^2 and differ only
   in coefficients, signs, and the final assembly), so they live in
   one file; the public interface keeps separate functions per
   mathematical function.

   Dispatch mirrors exp_rs.c: the pre32 variants use hardcoded
   straight-line routines (n <= 17 for atan/atanh, n <= 10 for the
   sin/cos families) and otherwise the generic fallbacks; the pre64
   variants use the windowed general routines when x < 2^-128 or n
   exceeds the hardcoded range (n <= 35 resp. n <= 22), and hardcoded
   routines in between.

   Error bounds contributing to the FIXED_*_RS_MAX_ERR macros:
   hardcoded routines <= 12 ulp; the windowed general routines report
   at most 12-13 for these series (as for exp, their full-window
   increments never fire: the first atan denominator block covers odd
   numbers up to 33, i.e. k <= 16, and the first factorial block
   k <= 10, in both cases past the point where the windows have
   shrunk); the fallbacks stay below 10.  Alternating-series bounds
   are two-sided, hyperbolic ones one-sided. */

#if FLINT_BITS == 64

#include "trig_rs_hard.inc"
#if FLINT_BITS == 64
#include "hand_mulhi.inc"
#endif

#define TMP_ALLOC_LIMBS(size) TMP_ALLOC((size) * sizeof(ulong))
#define PW(i) (xpow + ((i) - 1) * pslot)

static void
propagate_add(nn_ptr ptr, slong pos, slong top, ulong cy)
{
    while (cy != 0 && pos <= top)
    {
        ulong t = ptr[pos] + cy;
        cy = (t < cy);
        ptr[pos] = t;
        pos++;
    }
}

static void
propagate_sub(nn_ptr ptr, slong pos, slong top, ulong cy)
{
    while (cy != 0 && pos <= top)
    {
        ulong t = ptr[pos];
        ptr[pos] = t - cy;
        cy = (t < cy);
        pos++;
    }
}

static void
boundary_mul(nn_ptr sb, slong xn, slong wn, slong wn2,
    nn_srcptr pwm, slong pslot, slong lzm, nn_ptr tp, ulong * err)
{
    slong wb = xn - wn;
    slong wb2 = xn - wn2;
    slong n = wn + 1;
    slong gb, top, lo;

    flint_mpn_mulhigh_n(tp, sb + wb, pwm + pslot - n, n);

    gb = wb - lzm;
    top = xn - lzm;

    if (top < wb2)
    {
        flint_mpn_zero(sb + wb2, xn - wb2 + 1);
        return;
    }

    lo = FLINT_MAX(gb, wb2);
    flint_mpn_copyi(sb + lo, tp + (lo - gb), top - lo + 1);

    if (gb > wb2)
        flint_mpn_zero(sb + wb2, gb - wb2);
    if (top < xn)
        flint_mpn_zero(sb + top + 1, xn - top);

    if (wn == xn && lzm == 0)
        *err += 3;
}

/* write z = x^2 and its powers z^2..z^m into slots; power z^i gets its
   top G_i = xn + 1 - Lv*i limbs (one sub-grid guard limb) */
static void
write_z_powers(nn_ptr xpow, slong pslot, nn_srcptr x, slong xn,
    slong zL, slong Lv, slong m, nn_ptr tp)
{
    slong i, g, G, xsig = xn - zL;

    flint_mpn_zero(xpow, m * pslot);

    G = xn + 1 - Lv;
    g = FLINT_MIN(xsig, G);
    if (g >= 1)
    {
        flint_mpn_sqrhigh(tp, x + xsig - g, g);
        flint_mpn_copyi(PW(1) + pslot - g, tp, g);
    }

    for (i = 2; i <= m; i++)
    {
        G = xn + 1 - Lv * i;
        if (G < 1)
            continue;
        if (i % 2 == 0)
            flint_mpn_sqrhigh(PW(i) + pslot - G, PW(i / 2) + pslot - G, G);
        else
            flint_mpn_mulhigh_n(PW(i) + pslot - G,
                PW(i / 2) + pslot - G, PW(i / 2 + 1) + pslot - G, G);
    }
}

static void
_fixed_atan_rs_gen(nn_ptr y, ulong * error, nn_srcptr x, slong xn,
    ulong N, int alternating)
{
    slong zL, Lv, m, k, power, base, guard, nblocks, blk;
    slong wn, wb, wn2, drop, pslot, un;
    ulong c, new_denom, old_denom, cy, err, P, hi, lo;
    nn_ptr xpow, sbuf, sb, tp, xpad, pwp, bden;
    slong * bstart;
    TMP_INIT;

    zL = 0;
    while (zL < xn && x[xn - 1 - zL] == 0)
        zL++;

    if (zL == xn || N == 0)
    {
        flint_mpn_zero(y, xn);
        error[0] = 0;
        return;
    }

    FLINT_ASSERT(zL >= 1);
    Lv = 2 * zL;

    /* terms with 2k+1 >= (xn+1)/zL + 1 start below the sub-ulp guard */
    if (N > (ulong)(((xn + 1) / zL) / 2 + 1))
        N = (ulong)(((xn + 1) / zL) / 2 + 1);

    if (N == 1)
    {
        flint_mpn_copyi(y, x, xn);
        error[0] = 1;
        return;
    }

    m = 2;
    while (2 * m * m < (slong) N)
        m += 2;

    TMP_START;

    pslot = xn + 2;
    xpow = TMP_ALLOC_LIMBS(m * pslot + (xn + 4) + (xn + 4) + (xn + 1));
    sbuf = xpow + m * pslot;
    tp = sbuf + (xn + 4);
    xpad = tp + (xn + 4);
    bden = TMP_ALLOC_LIMBS(N + 2);
    bstart = (slong *) TMP_ALLOC((N + 2) * sizeof(slong));

    /* ascending prepass: blocks of consecutive odd numbers whose
       product fits a limb */
    nblocks = 0;
    P = 1;
    bstart[0] = 0;
    for (k = 1; k < (slong) N; k++)
    {
        umul_ppmm(hi, lo, P, (ulong)(2 * k + 1));
        if (hi != 0)
        {
            bden[nblocks++] = P;
            bstart[nblocks] = k;
            P = (ulong)(2 * k + 1);
        }
        else
            P = lo;
    }
    bden[nblocks++] = P;

    guard = (nblocks == 1) ? 1 : 2;

    flint_mpn_zero(sbuf, xn + 4);
    sb = sbuf + 1;

    write_z_powers(xpow, pslot, x, xn, zL, Lv, m, tp);

    power = (slong)(N - 1) % m;
    base = (slong)(N - 1) - power;
    drop = Lv * base - guard;
    wn = (drop <= 0) ? xn : FLINT_MAX(xn - drop, 1);
    wb = xn - wn;
    un = wn - Lv * power;
    pwp = xpow + power * pslot;

    blk = nblocks - 1;
    while (bstart[blk] > (slong) N - 1)
        blk--;

    err = 6 + 2 * (guard == 1);

    for (k = N - 1; k >= 0; k--)
    {
        if (k < bstart[blk])
        {
            old_denom = bden[blk];
            blk--;
            new_denom = bden[blk];

            if (alternating && (k % 2 == 0))
                sb[xn] += old_denom;

            sb[xn + 1] = mpn_mul_1(sb + wb, sb + wb, wn + 1, new_denom);
            mpn_divrem_1(sb + wb, 0, sb + wb, wn + 2, old_denom);

            if (alternating && (k % 2 == 0))
                sb[xn] -= new_denom;

            if (wn == xn)
                err += 2;
        }
        c = bden[blk] / (ulong)(2 * k + 1);

        if (power == 0)
        {
            if (alternating & k)
                sb[xn] -= c;
            else
                sb[xn] += c;

            if (k != 0)
            {
                base -= m;
                drop = Lv * base - guard;
                wn2 = (drop <= 0) ? xn : FLINT_MAX(xn - drop, 1);
                boundary_mul(sb, xn, wn, wn2, PW(m), pslot, Lv * m, tp, &err);
                wn = wn2;
                wb = xn - wn;
            }

            power = m - 1;
            un = wn - Lv * (m - 1);
            pwp = xpow + (m - 1) * pslot;
        }
        else
        {
            if (un > 0)
            {
                if (alternating & k)
                {
                    cy = mpn_submul_1(sb + wb, pwp - un, un, c);
                    if (un == wn)
                        sb[xn] -= cy;
                    else
                        propagate_sub(sb, wb + un, xn, cy);
                }
                else
                {
                    cy = mpn_addmul_1(sb + wb, pwp - un, un, c);
                    if (un == wn)
                        sb[xn] += cy;
                    else
                        propagate_add(sb, wb + un, xn, cy);
                }
            }
            un += Lv;
            pwp -= pslot;
            power--;
        }
    }

    mpn_divrem_1(sb, 0, sb, xn + 1, bden[0]);
    err += 2;

    xpad[0] = 0;
    flint_mpn_copyi(xpad + 1, x, xn);
    flint_mpn_mulhigh_n(tp, sb, xpad, xn + 1);
    flint_mpn_copyi(y, tp, xn);
    err += 2;

    error[0] = err;

    TMP_END;
}

static void
_fixed_sin_cos_rs_gen(nn_ptr ysin, nn_ptr ycos, ulong * error,
    nn_srcptr x, slong xn, ulong N, int sinonly, int alternating)
{
    slong zL, Lv, m, k, power, base, guard;
    slong wn, wb, wn2, drop, pslot, un;
    ulong c, dsf, F, cy, err, err1, hi, lo;
    nn_ptr xpow, sbuf, sb, tp, xpad, pwp;
    int cosorsin, single;
    TMP_INIT;

    zL = 0;
    while (zL < xn && x[xn - 1 - zL] == 0)
        zL++;

    if (zL == xn || N <= 1)
    {
        if (!sinonly)
        {
            flint_mpn_zero(ycos, xn);
            ycos[xn] = (N >= 1);
        }
        if (zL == xn || N == 0)
            flint_mpn_zero(ysin, xn);
        else
            flint_mpn_copyi(ysin, x, xn);
        ysin[xn] = 0;
        error[0] = 2;
        return;
    }

    FLINT_ASSERT(zL >= 1);
    Lv = 2 * zL;

    if (N > (ulong)(((xn + 1) / zL) / 2 + 1))
        N = (ulong)(((xn + 1) / zL) / 2 + 1);

    m = 2;
    while (2 * m * m < (slong) N)
        m += 2;

    /* single denominator block for both series iff (2N-1)! fits */
    single = (2 * N - 1 <= 20);
    guard = single ? 1 : 2;

    TMP_START;

    pslot = xn + 2;
    xpow = TMP_ALLOC_LIMBS(m * pslot + (xn + 4) + (xn + 4) + (xn + 1));
    sbuf = xpow + m * pslot;
    tp = sbuf + (xn + 4);
    xpad = tp + (xn + 4);

    write_z_powers(xpow, pslot, x, xn, zL, Lv, m, tp);

    err = 0;

    for (cosorsin = sinonly ? 1 : 0; cosorsin <= 1; cosorsin++)
    {
        flint_mpn_zero(sbuf, xn + 4);
        sb = sbuf + 1;

        power = (slong)(N - 1) % m;
        base = (slong)(N - 1) - power;
        drop = Lv * base - guard;
        wn = (drop <= 0) ? xn : FLINT_MAX(xn - drop, 1);
        wb = xn - wn;
        un = wn - Lv * power;
        pwp = xpow + power * pslot;

        err1 = 5 + 2 * (guard == 1);

        /* on-the-fly: c = f(N-1)!/f(k)!, dsf = f(N-1)!/f(k-1)!,
           f(k) = 2k + cosorsin; closing a block divides by dsf */
        c = 1;
        F = (2 * (N - 1) + cosorsin >= 2)
            ? (ulong)(2 * (N - 1) + cosorsin) * (2 * (N - 1) + cosorsin - 1)
            : 1;
        dsf = F;

        for (k = N - 1; k >= 0; k--)
        {
            if (power == 0)
            {
                if (alternating & k)
                    sb[xn] -= c;
                else
                    sb[xn] += c;

                if (k != 0)
                {
                    base -= m;
                    drop = Lv * base - guard;
                    wn2 = (drop <= 0) ? xn : FLINT_MAX(xn - drop, 1);
                    boundary_mul(sb, xn, wn, wn2, PW(m), pslot, Lv * m,
                        tp, &err1);
                    wn = wn2;
                    wb = xn - wn;
                }

                power = m - 1;
                un = wn - Lv * (m - 1);
                pwp = xpow + (m - 1) * pslot;
            }
            else
            {
                if (un > 0)
                {
                    if (alternating & k)
                    {
                        cy = mpn_submul_1(sb + wb, pwp - un, un, c);
                        if (un == wn)
                            sb[xn] -= cy;
                        else
                            propagate_sub(sb, wb + un, xn, cy);
                    }
                    else
                    {
                        cy = mpn_addmul_1(sb + wb, pwp - un, un, c);
                        if (un == wn)
                            sb[xn] += cy;
                        else
                            propagate_add(sb, wb + un, xn, cy);
                    }
                }
                un += Lv;
                pwp -= pslot;
                power--;
            }

            /* step coefficients k -> k-1: c_{k-1} = c_k * F(k) = dsf_k
               and dsf_{k-1} = dsf_k * F(k-1), F(j) = f(j)*(f(j)-1) */
            if (k > 0)
            {
                F = (2 * (k - 1) + cosorsin >= 2)
                    ? (ulong)(2 * (k - 1) + cosorsin)
                        * (2 * (k - 1) + cosorsin - 1)
                    : 1;
                umul_ppmm(hi, lo, dsf, F);
                if (hi != 0)
                {
                    /* close the block: divide window by dsf; in the
                       negative phase (last added term k odd) offset
                       the units limb to keep the dividend positive */
                    if (alternating && (k & 1))
                        sb[xn] += dsf;
                    mpn_divrem_1(sb + wb, 0, sb + wb, wn + 1, dsf);
                    if (alternating && (k & 1))
                        sb[xn] -= 1;
                    if (wn == xn)
                        err1 += 2;
                    c = 1;
                    dsf = F;
                }
                else
                {
                    c = dsf;
                    dsf = lo;
                }
            }
        }

        mpn_divrem_1(sb, 0, sb, xn + 1,
            dsf);
        err1 += 2;

        if (cosorsin == 0)
        {
            flint_mpn_copyi(ycos, sb, xn + 1);
        }
        else
        {
            xpad[0] = 0;
            flint_mpn_copyi(xpad + 1, x, xn);
            flint_mpn_mulhigh_n(tp, sb, xpad, xn + 1);
            flint_mpn_copyi(ysin, tp, xn + 1);
            err1 += 2;
        }

        err = FLINT_MAX(err, err1);
    }

    error[0] = err;

    TMP_END;
}

#endif /* FLINT_BITS == 64 */

/* leading zero bits of (x, n), or -1 if x = 0 */
static slong
_fixed_lzb(nn_srcptr x, slong n)
{
    slong i = n - 1;

    while (i >= 0 && x[i] == 0)
        i--;

    if (i < 0)
        return -1;

    return FLINT_BITS * (n - 1 - i) + flint_clz(x[i]);
}

/* sin/cos (alternating = 1) or sinh/cosh (alternating = 0) of (x, n)
   -> (ysin, n + 1) and (ycos, n + 1), either possibly NULL, for
   0 <= x < 2^-32, any n >= 1, at constant full precision with
   coefficients generated on the fly.  Descending numerators for the
   factorial series chain integrally: with f(k) = 2k + [odd series]
   and F(j) = f(j)(f(j)-1) (or 1 when f(j) < 2), the numerator for
   term k - 1 equals the running block divisor for term k.  For the
   alternating case, even m keeps the accumulator positive at every
   boundary multiplication, and block divisions in the negative phase
   (last added term of odd index) offset the units limb by the old
   divisor before dividing and subtract 1 after.  Error < 10 ulp. */
void
_fixed_sin_cos_rs_fallback(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int alternating)
{
    slong zb, N, m, j, k, power, cs;
    ulong c, dsf, F, cy;
    nn_ptr z, pw, s, tp, pad;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);

    zb = _fixed_lzb(x, n);

    /* the number of terms is chosen from the actual leading zero
       bits, so the wider range x < 2^-16 is served too (used by
       the wide-range family beyond its hardcoded reach) */
    if (zb >= 0)
        FLINT_ASSERT(zb >= 16);

    /* first dropped terms x^(2N)/(2N)!, x^(2N+1)/(2N+1)! below
       2^(-FLINT_BITS (n+1)) */
    N = (zb < 0) ? 1 : (FLINT_BITS * (n + 1)) / (2 * zb) + 1;

    if (N == 1)
    {
        if (ysin != NULL)
        {
            flint_mpn_copyi(ysin, x, n);
            ysin[n] = 0;
        }
        if (ycos != NULL)
        {
            flint_mpn_zero(ycos, n);
            ycos[n] = 1;
        }
        return;
    }

    m = 2;
    while (m * m < N)
        m += 2;

    TMP_START;
    z = TMP_ALLOC((n + ((m - 1) * n) + 3 * (n + 2)) * sizeof(ulong));
    pw = z + n;
    s = pw + (m - 1) * n;
    tp = s + (n + 2);
    pad = tp + (n + 2);

#define PWZ(i) ((i) == 1 ? (nn_srcptr) z : (nn_srcptr) (pw + ((i) - 2) * n))

    flint_mpn_sqrhigh(z, x, n);

    for (j = 2; j <= FLINT_MIN(m, N - 1); j++)
    {
        if (j % 2 == 0)
            flint_mpn_sqrhigh(pw + (j - 2) * n, PWZ(j / 2), n);
        else
            flint_mpn_mulhigh_n(pw + (j - 2) * n, PWZ(j / 2),
                PWZ(j / 2 + 1), n);
    }

    for (cs = 0; cs <= 1; cs++)
    {
        nn_ptr out = (cs == 0) ? ycos : ysin;

        if (out == NULL)
            continue;

        flint_mpn_zero(s, n + 2);
        power = (N - 1) % m;

        c = 1;
        F = (2 * (N - 1) + cs >= 2)
            ? (ulong)(2 * (N - 1) + cs) * (2 * (N - 1) + cs - 1) : 1;
        dsf = F;

        for (k = N - 1; k >= 0; k--)
        {
            if (power == 0)
            {
                if (alternating & k)
                    s[n] -= c;
                else
                    s[n] += c;

                if (k != 0)
                {
                    /* multiply the full accumulator by z^m */
                    pad[0] = 0;
                    flint_mpn_copyi(pad + 1, PWZ(m), n);
                    flint_mpn_mulhigh_n(tp, s, pad, n + 1);
                    flint_mpn_copyi(s, tp, n + 1);
                }

                power = m - 1;
            }
            else
            {
                if (alternating & k)
                {
                    cy = mpn_submul_1(s, PWZ(power), n, c);
                    s[n] -= cy;
                }
                else
                {
                    cy = mpn_addmul_1(s, PWZ(power), n, c);
                    s[n] += cy;
                }
                power--;
            }

            /* step coefficients k -> k-1: c_{k-1} = dsf_k and
               dsf_{k-1} = dsf_k F(k-1); close the block by dividing
               when the divisor would overflow */
            if (k > 0)
            {
                ulong hi, lo;

                F = (2 * (k - 1) + cs >= 2)
                    ? (ulong)(2 * (k - 1) + cs) * (2 * (k - 1) + cs - 1)
                    : 1;
                umul_ppmm(hi, lo, dsf, F);
                if (hi != 0)
                {
                    if (alternating && (k & 1))
                        s[n] += dsf;
                    mpn_divrem_1(s, 0, s, n + 1, dsf);
                    if (alternating && (k & 1))
                        s[n] -= 1;
                    c = 1;
                    dsf = F;
                }
                else
                {
                    c = dsf;
                    dsf = lo;
                }
            }
        }

        mpn_divrem_1(s, 0, s, n + 1, dsf);

        if (cs == 0)
            flint_mpn_copyi(out, s, n + 1);
        else
        {
            /* multiply by x */
            pad[0] = 0;
            flint_mpn_copyi(pad + 1, x, n);
            flint_mpn_mulhigh_n(tp, s, pad, n + 1);
            flint_mpn_copyi(out, tp, n + 1);
        }
    }

#undef PWZ

    TMP_END;
}

/* atan (alternating = 1) or atanh (alternating = 0) of (x, n)
   -> (res, n), for 0 <= x < 2^-32, any n >= 1, at constant full
   precision.  Odd reciprocals do not chain integrally downward, so an
   ascending prepass collects blocks of consecutive odd numbers whose
   product fits a limb; the descending evaluation then uses one exact
   division D/(2k+1) per term and the classical rescale (multiply by
   the next block divisor, divide by the previous, with units-limb
   offsets in the negative phase) at block changes.  Error < 10 ulp. */
void
_fixed_atan_rs_fallback(nn_ptr res, nn_srcptr x, slong n, int alternating)
{
    slong zb, N, m, j, k, power, blk, nblocks;
    ulong c, dsf, cy, new_denom, old_denom;
    nn_ptr z, pw, s, tp, pad, bden;
    slong * bstart;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);

    zb = _fixed_lzb(x, n);

    if (zb >= 0)
        FLINT_ASSERT(zb >= 16);

    N = (zb < 0) ? 1 : (FLINT_BITS * (n + 1)) / (2 * zb) + 1;

    if (N == 1)
    {
        flint_mpn_copyi(res, x, n);
        return;
    }

    m = 2;
    while (m * m < N)
        m += 2;

    TMP_START;
    z = TMP_ALLOC((n + ((m - 1) * n) + 3 * (n + 2)) * sizeof(ulong));
    pw = z + n;
    s = pw + (m - 1) * n;
    tp = s + (n + 2);
    pad = tp + (n + 2);
    bden = TMP_ALLOC((N + 2) * sizeof(ulong));
    bstart = (slong *) TMP_ALLOC((N + 2) * sizeof(slong));

#define PWZ(i) ((i) == 1 ? (nn_srcptr) z : (nn_srcptr) (pw + ((i) - 2) * n))

    flint_mpn_sqrhigh(z, x, n);

    for (j = 2; j <= FLINT_MIN(m, N - 1); j++)
    {
        if (j % 2 == 0)
            flint_mpn_sqrhigh(pw + (j - 2) * n, PWZ(j / 2), n);
        else
            flint_mpn_mulhigh_n(pw + (j - 2) * n, PWZ(j / 2),
                PWZ(j / 2 + 1), n);
    }

    /* ascending prepass: blocks of consecutive odd numbers whose
       product fits a limb */
    nblocks = 0;
    dsf = 1;
    bstart[0] = 0;
    for (k = 1; k < N; k++)
    {
        ulong hi, lo;

        umul_ppmm(hi, lo, dsf, (ulong)(2 * k + 1));
        if (hi != 0)
        {
            bden[nblocks++] = dsf;
            bstart[nblocks] = k;
            dsf = (ulong)(2 * k + 1);
        }
        else
            dsf = lo;
    }
    bden[nblocks++] = dsf;

    flint_mpn_zero(s, n + 3);
    power = (N - 1) % m;

    blk = nblocks - 1;

    for (k = N - 1; k >= 0; k--)
    {
        if (k < bstart[blk])
        {
            old_denom = bden[blk];
            blk--;
            new_denom = bden[blk];

            if (alternating && (k % 2 == 0))
                s[n] += old_denom;

            s[n + 1] = mpn_mul_1(s, s, n + 1, new_denom);
            mpn_divrem_1(s, 0, s, n + 2, old_denom);

            if (alternating && (k % 2 == 0))
                s[n] -= new_denom;
        }

        c = bden[blk] / (ulong)(2 * k + 1);

        if (power == 0)
        {
            if (alternating & k)
                s[n] -= c;
            else
                s[n] += c;

            if (k != 0)
            {
                pad[0] = 0;
                flint_mpn_copyi(pad + 1, PWZ(m), n);
                flint_mpn_mulhigh_n(tp, s, pad, n + 1);
                flint_mpn_copyi(s, tp, n + 1);
            }

            power = m - 1;
        }
        else
        {
            if (alternating & k)
            {
                cy = mpn_submul_1(s, PWZ(power), n, c);
                s[n] -= cy;
            }
            else
            {
                cy = mpn_addmul_1(s, PWZ(power), n, c);
                s[n] += cy;
            }
            power--;
        }
    }

    mpn_divrem_1(s, 0, s, n + 1, bden[0]);

    pad[0] = 0;
    flint_mpn_copyi(pad + 1, x, n);
    flint_mpn_mulhigh_n(tp, s, pad, n + 1);
    flint_mpn_copyi(res, tp, n);

#undef PWZ

    TMP_END;
}

/* dispatch helpers ********************************************************/

static void
_fixed_sc_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n,
    int alternating)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT((x[n - 1] >> (FLINT_BITS - 32)) == 0);

#if FLINT_BITS == 64
    if (x[n - 1] != 0)
    {
        /* 2^-64 <= x < 2^-32 */
        if (n <= 10)
        {
            (alternating ? _fixed_sin_cos_rs32_tab
                         : _fixed_sinh_cosh_rs32_tab)[n](ysin, ycos, x);
            return;
        }
    }
    else
    {
        ulong error;

        if ((n >= 2 && x[n - 2] == 0) || n > 22)
        {
            if (ysin == NULL)
            {
                /* the general routine has no cos-only mode */
                nn_ptr tsin;
                TMP_INIT;
                TMP_START;
                tsin = TMP_ALLOC((n + 1) * sizeof(ulong));
                _fixed_sin_cos_rs_gen(tsin, ycos, &error, x, n,
                    (ulong) n + 2, 0, alternating);
                TMP_END;
            }
            else
                _fixed_sin_cos_rs_gen(ysin, ycos, &error, x, n,
                    (ulong) n + 2, ycos == NULL, alternating);
        }
        else
            (alternating ? _fixed_sin_cos_rs_tab
                         : _fixed_sinh_cosh_rs_tab)[n](ysin, ycos, x);
        return;
    }
#endif
    _fixed_sin_cos_rs_fallback(ysin, ycos, x, n, alternating);
}

static void
_fixed_at_rs(nn_ptr res, nn_srcptr x, slong n, int alternating)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT((x[n - 1] >> (FLINT_BITS - 32)) == 0);

#if FLINT_BITS == 64
    if (x[n - 1] != 0)
    {
        /* 2^-64 <= x < 2^-32 */
        if (n <= 17)
        {
            (alternating ? _fixed_atan_rs32_tab
                         : _fixed_atanh_rs32_tab)[n](res, x);
            return;
        }
    }
    else
    {
        ulong error;

        if ((n >= 2 && x[n - 2] == 0) || n > 35)
            _fixed_atan_rs_gen(res, &error, x, n, (ulong) n + 2,
                alternating);
        else
            (alternating ? _fixed_atan_rs_tab
                         : _fixed_atanh_rs_tab)[n](res, x);
        return;
    }
#endif
    _fixed_atan_rs_fallback(res, x, n, alternating);
}

/* public functions ********************************************************/

void fixed_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n)
{ _fixed_sc_rs(ysin, ycos, x, n, 1); }

void fixed_sin_rs(nn_ptr res, nn_srcptr x, slong n)
{ _fixed_sc_rs(res, NULL, x, n, 1); }

void fixed_cos_rs(nn_ptr res, nn_srcptr x, slong n)
{ _fixed_sc_rs(NULL, res, x, n, 1); }

void fixed_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, nn_srcptr x, slong n)
{ _fixed_sc_rs(ysinh, ycosh, x, n, 0); }

void fixed_sinh_rs(nn_ptr res, nn_srcptr x, slong n)
{ _fixed_sc_rs(res, NULL, x, n, 0); }

void fixed_cosh_rs(nn_ptr res, nn_srcptr x, slong n)
{ _fixed_sc_rs(NULL, res, x, n, 0); }

void fixed_atan_rs(nn_ptr res, nn_srcptr x, slong n)
{ _fixed_at_rs(res, x, n, 1); }

void fixed_atanh_rs(nn_ptr res, nn_srcptr x, slong n)
{ _fixed_at_rs(res, x, n, 0); }
