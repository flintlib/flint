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

/* Rectangular splitting evaluation of exp on [0, 2^-32) resp.
   [0, 2^-64), following "Fast multiple precision exp(x)
   with precomputations" (van der Hoeven and Johansson) and "Efficient
   implementation of elementary functions in the medium-precision
   range" (Johansson).

   Dispatch: fixed_exp_rs_pre32 uses a hardcoded straight-line routine
   for n <= 11 on 64-bit machines (the hardcoded table is indexed by
   the number of Taylor terms N with n = ceil(N/2); N = 2n keeps the
   dropped tail below one sub-ulp, and N = 21 serves n = 11 since
   21! > 2^33 damps the half-limb overhang), and otherwise the generic
   fallback, since the windowed general routine requires a whole
   leading zero limb.  fixed_exp_rs_pre64 checks for a second leading
   zero limb (x < 2^-128) and then calls the windowed general routine
   _fixed_exp_rs_gen, which exploits the extra zeros; otherwise a
   hardcoded routine for n <= 21, and the general routine beyond.

   Error bounds contributing to FIXED_EXP_RS_MAX_ERR: hardcoded
   routines <= 8 ulp plus a sub-ulp tail; _fixed_exp_rs_gen reports at
   most 7 for exp (its +2 full-window division and +3 full-window
   boundary increments never fire, since the first coefficient block
   close happens at k = 20 where the window has already shrunk, and
   zL * m >= 1); the generic fallback is bounded below.  All
   truncations are downward, so the bound is one-sided. */

#if FLINT_BITS == 64

#include "exp_rs_hard.inc"

#define TMP_ALLOC_LIMBS(size) TMP_ALLOC((size) * sizeof(ulong))
#define GUARD_LIMBS 2
#define PW(i) (xpow + ((i) - 1) * pslot)   /* slot of x^i, 1-indexed */

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

/* multiply the current window [wb, xn] by x^m (top limbs of its slot,
   zero-extended one limb below when the read overruns), writing the
   result aligned to the new window [wb2, xn] */
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

static void
_fixed_exp_rs_gen(nn_ptr y, ulong * error, nn_srcptr x, slong xn, ulong N)
{
    slong zL, m, i, k, power, base, guard;
    slong wn, wb, wn2, drop, pslot, un, G;
    ulong c, dsf, cy, err;
    nn_ptr xpow, sbuf, sb, tp, pwp;
    TMP_INIT;

    /* leading zero limbs of x; the contract requires x < 2^-64 */
    zL = 0;
    while (zL < xn && x[xn - 1 - zL] == 0)
        zL++;

    if (zL == xn || N <= 1)             /* x = 0 or y = 1 (or empty sum) */
    {
        flint_mpn_zero(y, xn);
        y[xn] = (N >= 1);
        error[0] = 0;
        return;
    }

    FLINT_ASSERT(zL >= 1);

    /* drop terms whose first limb is below the guard limb under the
       final ulp: x^k < 2^(-64*zL*k) <= 2^(-64*(xn+1)) for
       k >= (xn+1)/zL; the dropped tail is < 2^(-63*...) << 1 ulp */
    if (N > (ulong)(xn + 1) / zL + 1)
        N = (ulong)(xn + 1) / zL + 1;

    if (N == 2)
    {
        flint_mpn_copyi(y, x, xn);
        y[xn] = 1;
        error[0] = 1;
        return;
    }

    /* rectangular splitting parameter: even m ~ sqrt(N/2); the power
       computations cost full-length mulhighs while the boundary
       multiplications shrink with the window, so the optimum sits
       below sqrt(N) (measured within a few percent of per-size optima
       across xn = 10..320, zL = 1..3) */
    m = 2;
    while (2 * m * m < (slong) N)
        m += 2;

    /* Window guard limbs.  Every truncation event committed at a
       shrunk window (term reads, boundary mulhighs, divisions) is
       later multiplied by the pending x^base <= 2^(-64*zL*base) while
       the window bottom sits zL*base - guard limbs below the top, so
       each is suppressed by 2^(-64*guard) regardless of the block
       divisors.  With N <= 21 there is a single coefficient block
       (divisor (N-1)!) and hence only the final division; guard = 1
       then gives the budget: full-window term drops sum_{k} 1/k!
       <= 1.72, last boundary <= 3/(N-1)! <= 0.13, final division
       <= 1, total < 4 ulp. */
    guard = (N <= 21) ? 1 : GUARD_LIMBS;
    /* keep at least one limb of x^m on the extended grid; N >= 3
       already guarantees this for m = 2 */
    while (m > 2 && xn + 1 - zL * m < 1)
        m -= 2;

    TMP_START;

    pslot = xn + 2;
    xpow = TMP_ALLOC_LIMBS(m * pslot + (xn + 3) + (xn + 4));
    sbuf = xpow + m * pslot;
    tp = sbuf + (xn + 3);

    /* power slots need no zeroing: every read stays within the limbs
       written below, except the boundary multiplication's one-limb
       zero-extension below x^m's stored region */
    flint_mpn_zero(sbuf, xn + 3);
    sb = sbuf + 1;
    /* the boundary reads extend up to guard limbs below x^m's stored
       region (n = wn + 1 <= G_m + guard); zero exactly those */
    flint_mpn_zero(PW(m) + pslot - (xn + 1 - zL * m) - guard, guard);

    /* x, then x^2..x^m in place (all shift deltas are zero) */
    flint_mpn_copyi(PW(1) + pslot - (xn - zL), x, xn - zL);
    for (i = 2; i <= m; i++)
    {
        G = xn + 1 - zL * i;
        if (G < 1)
            continue;
        if (i % 2 == 0)
            flint_mpn_sqrhigh(PW(i) + pslot - G,
                PW(i / 2) + pslot - G, G);
        else
            flint_mpn_mulhigh_n(PW(i) + pslot - G,
                PW(i / 2) + pslot - G, PW(i / 2 + 1) + pslot - G, G);
    }

    power = (slong)(N - 1) % m;
    base = (slong)(N - 1) - power;
    drop = zL * base - guard;
    wn = (drop <= 0) ? xn : FLINT_MAX(xn - drop, 1);
    wb = xn - wn;
    un = wn - zL * power;
    pwp = xpow + power * pslot;

    err = (guard == 1) ? 2 : 5;

    /* running numerator and block divisor */
    c = 1;
    dsf = N - 1;

    for (k = N - 1; k >= 0; k--)
    {
        if (k < (slong) N - 1)
        {
            ulong kk = (k >= 1) ? (ulong) k : 1;
            ulong hi, lo;

            umul_ppmm(hi, lo, dsf, kk);
            if (hi != 0)
            {
                /* close the block [k+1, ...]: divide the window */
                mpn_divrem_1(sb + wb, 0, sb + wb, wn + 1, dsf);
                if (wn == xn)
                    err += 2;
                c = 1;
                dsf = kk;
            }
            else
            {
                c = c * (ulong)(k + 1);
                dsf = lo;
            }
        }

        if (power == 0)
        {
            sb[xn] += c;

            if (k != 0)
            {
                base -= m;
                drop = zL * base - guard;
                wn2 = (drop <= 0) ? xn : FLINT_MAX(xn - drop, 1);
                boundary_mul(sb, xn, wn, wn2, PW(m), pslot, zL * m, tp, &err);
                wn = wn2;
                wb = xn - wn;
            }

            power = m - 1;
            un = wn - zL * (m - 1);
            pwp = xpow + (m - 1) * pslot;
        }
        else
        {
            if (un > 0)
            {
                cy = mpn_addmul_1(sb + wb, pwp - un, un, c);
                if (un == wn)
                    sb[xn] += cy;
                else
                    propagate_add(sb, wb + un, xn, cy);
            }
            un += zL;
            pwp -= pslot;
            power--;
        }
    }

    mpn_divrem_1(sb, 0, sb, xn + 1, dsf);
    err += 2;

    flint_mpn_copyi(y, sb, xn + 1);
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

/* exp((x, n)) -> (res, n + 1) for 0 <= x < 2^-32, any n >= 1, at
   constant full precision with coefficients generated on the fly (as
   in _arb_exp_taylor_rs, but with no table limit on N).  The number
   of terms is chosen from the actual leading zero bits of x.

   Error (one-sided; every truncation is downward): the powers
   x^2..x^m accumulate a few ulp each through the mulhigh chain but
   enter the sum with weight at most 1/k!, contributing below 2 ulp in
   total; each boundary multiplication truncates at most 2 ulp, damped
   by the pending power of x except for the last one; block divisions
   floor at most 1 ulp each, all but the last damped by the following
   divisors.  Total < 6 ulp. */
void
_fixed_exp_rs_fallback(nn_ptr res, nn_srcptr x, slong n)
{
    slong zb, N, m, j, k, power;
    ulong c, dsf, cy;
    nn_ptr pw, s, tp, xmpad;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);

    zb = _fixed_lzb(x, n);

    if (zb < 0)         /* x = 0 */
    {
        flint_mpn_zero(res, n);
        res[n] = 1;
        return;
    }

    FLINT_ASSERT(zb >= 32);

    /* first dropped term x^N/N! < 2^(-FLINT_BITS (n+1)) */
    N = (FLINT_BITS * (n + 1) + zb - 1) / zb;

    m = 2;
    while (m * m < N)
        m++;
    m = FLINT_MIN(m, N - 1);

    TMP_START;
    pw = TMP_ALLOC((((m - 1) * n) + 3 * (n + 2)) * sizeof(ulong));
    s = pw + (m - 1) * n;
    tp = s + (n + 2);
    xmpad = tp + (n + 2);

#define PWX(i) ((i) == 1 ? x : (nn_srcptr) (pw + ((i) - 2) * n))

    for (j = 2; j <= m; j++)
    {
        if (j % 2 == 0)
            flint_mpn_sqrhigh(pw + (j - 2) * n, PWX(j / 2), n);
        else
            flint_mpn_mulhigh_n(pw + (j - 2) * n, PWX(j / 2),
                PWX(j / 2 + 1), n);
    }

    flint_mpn_zero(s, n + 2);

    power = (N - 1) % m;

    /* running numerator and block divisor */
    c = 1;
    dsf = N - 1;

    for (k = N - 1; k >= 0; k--)
    {
        if (k < N - 1)
        {
            ulong kk = (k >= 1) ? (ulong) k : 1;
            ulong hi, lo;

            umul_ppmm(hi, lo, dsf, kk);
            if (hi != 0)
            {
                /* close the coefficient block: divide by dsf */
                mpn_divrem_1(s, 0, s, n + 1, dsf);
                c = 1;
                dsf = kk;
            }
            else
            {
                c = c * (ulong)(k + 1);
                dsf = lo;
            }
        }

        if (power == 0)
        {
            s[n] += c;

            if (k != 0)
            {
                /* multiply the full accumulator by x^m */
                xmpad[0] = 0;
                flint_mpn_copyi(xmpad + 1, PWX(m), n);
                flint_mpn_mulhigh_n(tp, s, xmpad, n + 1);
                flint_mpn_copyi(s, tp, n + 1);
            }

            power = m - 1;
        }
        else
        {
            cy = mpn_addmul_1(s, PWX(power), n, c);
            s[n] += cy;
            power--;
        }
    }

    mpn_divrem_1(s, 0, s, n + 1, dsf);
    flint_mpn_copyi(res, s, n + 1);

#undef PWX

    TMP_END;
}

void
fixed_exp_rs(nn_ptr res, nn_srcptr x, slong n)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT((x[n - 1] >> (FLINT_BITS - 32)) == 0);

#if FLINT_BITS == 64
    if (x[n - 1] != 0)
    {
        /* 2^-64 <= x < 2^-32 */
        if (n <= 11)
        {
            _fixed_exp_rs32_tab[FLINT_MIN(2 * n, 21)](res, x);
            return;
        }
    }
    else
    {
        ulong error;

        if (n >= 2 && x[n - 2] == 0)
        {
            /* x < 2^-128: the windowed general routine exploits the
               extra leading zero limbs */
            _fixed_exp_rs_gen(res, &error, x, n, (ulong) n + 1);
        }
        else if (n <= 21)
            _fixed_exp_rs_tab[n](res, x);
        else
            _fixed_exp_rs_gen(res, &error, x, n, (ulong) n + 1);
        return;
    }
#endif
    _fixed_exp_rs_fallback(res, x, n);
}
