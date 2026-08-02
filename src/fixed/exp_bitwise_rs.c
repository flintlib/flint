/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "arb.h"
#include "fixed.h"

/* exp via bitwise argument reduction: subtract in turn each
   L_i = log(1 + 2^-i), i = 0, 1, ..., r, for which L_i <= x, evaluate
   the Taylor series on the reduced argument below 2^-r, and multiply
   back the used (1 + 2^-i) factors, each of which is a single
   shift-and-add.

   The classical single-pass reduction is sound on [0, 1): since
   (1 + 2^-i)^2 > 1 + 2^(1-i) we have L_{i-1} < 2 L_i, so by induction
   the remainder is below L_i after step i (base case 1 < 2 log 2),
   and each L_i is used at most once.  The stored logarithms are
   truncations, which lets the remainder creep at most one guard ulp
   per step above the exact invariant; a single extra conditional
   subtraction of L_r after the loop therefore guarantees a remainder
   strictly below L_r < 2^-r.

   All work happens at the output precision of n limbs -- callers
   needing sub-ulp accuracy are expected to pad the precision
   themselves, so padding internally would double up.  The error
   therefore grows linearly with r: at most r + 2 table entries are
   subtracted, each a truncation short of the true logarithm by less
   than two ulp (multiplying the result by exp of the deficit), and
   each of the at most r + 2 reconstruction shifts truncates below
   one ulp, all amplified by the remaining product < e; with the
   series error this gives less than 8.3 (r + 2) + 80 ulp, bounded by
   FIXED_EXP_BITWISE_RS_MAX_ERR(n, r) = 9 r + 100.  (With one guard
   limb these terms all drop below one output ulp -- measured 1.00
   ulp flat for r up to 320 -- so a caller-side padded limb absorbs
   the bound with 50+ bits to spare.)  r is clamped to
   FLINT_BITS n - 16: this keeps the reduction slop (r + 2 ulp)
   below L_r, so the single extra subtraction suffices, and past
   that point the series would degenerate anyway.

   The log table is generated on demand in two tiers -- binary
   splitting through arb for small i, where few series terms are
   shared between entries, and a fixed-point mpn multi-summation for
   large i, where one reciprocal per odd k serves EVERY entry (see
   _fixed_exp_logs_ensure below) -- and cached per thread at the
   largest precision and index range requested so far.  Each cached
   entry carries one guard limb below the value limbs; smaller
   requests read the top limbs of the cached entries, which equal
   the truncations at the smaller precision exactly. */

FLINT_TLS_PREFIX nn_ptr _fixed_exp_logs = NULL;
FLINT_TLS_PREFIX slong _fixed_exp_logs_n = 0;   /* limbs per entry,
                                                   INCLUDING one guard
                                                   limb at the bottom */
FLINT_TLS_PREFIX slong _fixed_exp_logs_r = 0;   /* entries 0..r */
static FLINT_TLS_PREFIX int _fixed_exp_logs_cleanup_registered = 0;

static void
_fixed_exp_logs_cleanup(void)
{
    flint_free(_fixed_exp_logs);
    _fixed_exp_logs = NULL;
    _fixed_exp_logs_n = 0;
    _fixed_exp_logs_r = 0;
    _fixed_exp_logs_cleanup_registered = 0;
}

#define LOGS_L(i) (_fixed_exp_logs + (i) * _fixed_exp_logs_n)
#define TAB_E(i) (tab + (i) * nc)


/* acc (nc limbs) += or -= t >> s, where t has nc limbs and the bits
   shifted below position 0 are dropped (they fall below the guard
   limb).  scratch holds nc limbs. */
static void
_fixed_tab_addshift(nn_ptr acc, nn_srcptr t, slong nc, slong s,
    int sub, nn_ptr scratch)
{
    slong q = s / FLINT_BITS, m = nc - q;
    int b = (int) (s % FLINT_BITS);
    nn_srcptr v;

    if (m <= 0)
        return;

    if (b == 0)
        v = t + q;
    else
    {
        mpn_rshift(scratch, t + q, m, b);
        v = scratch;
    }

    if (sub)
    {
        /* subtracted terms round UP (one extra guard ulp stands in
           for the dropped bits): together with the blanket bias
           applied to the leading term this keeps every entry AT OR
           BELOW the exact value, preserving the strict tab_i < 2^-i
           invariant of the reductions */
        mpn_sub(acc, acc, nc, v, m);
        mpn_sub_1(acc, acc, nc, 1);
    }
    else
        mpn_add(acc, acc, nc, v, m);
}

/* acc (nc limbs) -= 2^e, e < FLINT_BITS nc */
static void
_fixed_tab_subbit(nn_ptr acc, slong nc, slong e)
{
    slong q = e / FLINT_BITS;
    mpn_sub_1(acc + q, acc + q, nc - q,
        UWORD(1) << (e % FLINT_BITS));
}

/* The table entries L_i = log(1 + 2^-i), each nc limbs (the value
   scaled by 2^(FLINT_BITS nc), whose bottom limb is a guard).  Two
   tiers:

   Small i (below about sqrt(precision) bits, where the series is
   long): binary splitting.  For i <= 6 the arguments (2^i + 1)/2^i
   factor over the primes up to 17, so seven bsplit logarithms serve
   all of them; beyond that, fixed_log1p_2mexp_ui_bs (tab_bsplit.c)
   splits the series natively in mpn arithmetic, writing straight
   into the entry.  TODO: the remaining arb dependency of this tier
   is the seven-logarithm prime-vector combination for i <= 6.

   Large i: fixed-point multi-summation of

       log(1 + 2^-i) = sum_{k >= 1} (-1)^(k+1) 2^(-ik) / k

   grouped by the odd part of k: for odd k the even-index terms
   k' = k 2^l satisfy 2^(-ik')/k' = 2^(-(i k 2^l + l)) / k, so ONE
   reciprocal t = floor(2^wp / k) serves the whole group -- and,
   crucially, serves it for EVERY i at once, so the number of
   divisions is about wp / (2 iter_start) for the entire tier rather
   than per entry.  Everything else is shifted additions of t.  Each
   floored term errs below one guard ulp; a few thousand terms leave
   the value limbs equal to the exact truncation of L_i except when
   L_i lies within about 2^-50 of a limb boundary, and then by one
   ulp at most, which the reduction's accounting absorbs on either
   side. */
void
_fixed_exp_logs_ensure(slong nv, slong rc)
{
    slong nc, i, k, l, prec, iter_start;
    nn_ptr t, scratch, num;

    if (nv + 1 <= _fixed_exp_logs_n && rc <= _fixed_exp_logs_r)
        return;

    nc = FLINT_MAX(nv + 1, _fixed_exp_logs_n);
    rc = FLINT_MAX(rc, _fixed_exp_logs_r);

    flint_free(_fixed_exp_logs);
    _fixed_exp_logs = flint_malloc((rc + 1) * nc * sizeof(ulong));
    _fixed_exp_logs_n = nc;
    _fixed_exp_logs_r = rc;

    prec = FLINT_BITS * (nc - 1);

    if (prec <= 65536)
        iter_start = FLINT_MIN(n_sqrt(prec), rc + 1);
    else
        iter_start = rc + 1;

    /* tier 1: binary splitting */
    {
        arb_ptr lp;
        arb_t x;
        fmpz_t p, q;
        slong wp = FLINT_BITS * nc + 30;

        arb_init(x);
        fmpz_init(p);
        fmpz_init(q);
        lp = _arb_vec_init(7);

        arb_log_primes_vec_bsplit(lp, 7, wp);

        for (i = 0; i < iter_start && i <= rc; i++)
        {
            switch (i)
            {
                case 0:                 /* log 2 */
                    arb_set(x, lp + 0);
                    break;
                case 1:                 /* log(3/2) */
                    arb_sub(x, lp + 1, lp + 0, wp);
                    break;
                case 2:                 /* log(5/4) */
                    arb_mul_2exp_si(x, lp + 0, 1);
                    arb_sub(x, lp + 2, x, wp);
                    break;
                case 3:                 /* log(9/8) */
                    arb_mul_2exp_si(x, lp + 1, 1);
                    arb_submul_ui(x, lp + 0, 3, wp);
                    break;
                case 4:                 /* log(17/16) */
                    arb_mul_2exp_si(x, lp + 0, 2);
                    arb_sub(x, lp + 6, x, wp);
                    break;
                case 5:                 /* log(33/32) = log 3 + log 11
                                           - 5 log 2 */
                    arb_add(x, lp + 1, lp + 4, wp);
                    arb_submul_ui(x, lp + 0, 5, wp);
                    break;
                case 6:                 /* log(65/64) = log 5 + log 13
                                           - 6 log 2 */
                    arb_add(x, lp + 2, lp + 5, wp);
                    arb_submul_ui(x, lp + 0, 6, wp);
                    break;
                default:
                    /* native mpn binary splitting, straight into the
                       entry (atanh form, switching to the direct
                       series at high precision; see tab_bsplit.c) */
                    fixed_log1p_2mexp_ui_bs(LOGS_L(i), (ulong) i, nc);
                    continue;
            }

            if (!_fixed_tab_store_floor(LOGS_L(i), x, nc, wp))
                _fixed_tab_entry_exact(LOGS_L(i), 0, (ulong) i, nc);
        }

        _arb_vec_clear(lp, 7);
        arb_clear(x);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    /* tier 2: fixed-point multi-summation, directly in the target
       storage */
    if (iter_start <= rc)
    {
        slong wp = FLINT_BITS * nc;

        t = flint_malloc((nc + 2) * sizeof(ulong));
        scratch = flint_malloc(nc * sizeof(ulong));
        num = flint_malloc((nc + 1) * sizeof(ulong));

        /* the k = 1 group is pure bits: the leading term 2^-i and
           the corrections -2^(-(i 2^l + l)) for its even multiples */
        for (i = iter_start; i <= rc; i++)
        {
            nn_ptr acc = LOGS_L(i);

            flint_mpn_zero(acc, nc);
            acc[(wp - i) / FLINT_BITS] =
                UWORD(1) << ((wp - i) % FLINT_BITS);
            /* blanket bias: the series terms too small to represent
               (below the guard limb) sum to less than two guard
               ulps; over-subtracting them keeps the entry one-sided */
            mpn_sub_1(acc, acc, nc, 2);

            for (l = 1; ((i) << l) + l < wp; l++)
                _fixed_tab_subbit(acc, nc, wp - (((slong) i << l) + l));
        }

        /* odd k >= 3: one reciprocal per k serves every i */
        for (k = 3; k * iter_start < wp; k += 2)
        {
            flint_mpn_zero(num, nc);
            num[nc] = 1;                    /* 2^wp */
            mpn_divrem_1(t, 0, num, nc + 1, (ulong) k);
            FLINT_ASSERT(t[nc] == 0);

            for (i = iter_start; i <= rc && i * k < wp; i++)
            {
                nn_ptr acc = LOGS_L(i);

                _fixed_tab_addshift(acc, t, nc, i * k, 0, scratch);

                for (l = 1; ((i * k) << l) + l < wp; l++)
                    _fixed_tab_addshift(acc, t, nc,
                        (((slong) i * k) << l) + l, 1, scratch);
            }
        }

        flint_free(t);
        flint_free(scratch);
        flint_free(num);
    }

    if (!_fixed_exp_logs_cleanup_registered)
    {
        flint_register_cleanup_function(_fixed_exp_logs_cleanup);
        _fixed_exp_logs_cleanup_registered = 1;
    }
    /* Direct tail band: for 3i >= 64 nc every term of
       log(1 + 2^-i) = 2^-i - 2^(-2i-1) + 2^(-3i)/3 - ... from the
       third on falls below the entry, with the omitted alternating
       mass s satisfying 0 < s 2^(64 nc) <= 1/3, so the exact floor
       is a pure bit pattern: 2^(64 nc - i) - 2^(64 nc - 2i - 1)
       when the second term is still representable, and
       2^(64 nc - i) - 1 when it too falls below (then the omitted
       mass sits in (0, 1/2]).  Writing these directly covers the
       entries whose guard limb legitimately sits near wrapping.

       Exactness rescan below the band: the fast tiers produce
       entries at most a few guard ulps below the truth, so a guard
       limb within FIXED_TAB_GUARD_SLACK of wrapping means that
       deficit may have borrowed into the value limbs -- recompute
       such (now genuinely astronomically rare) entries via arb.
       Together with the exact-floor tier-1 store, every entry's
       value limbs equal the exact floor at their precision, hence
       at every shorter truncation as well. */
    {
        slong wp = FLINT_BITS * nc, band = (wp + 2) / 3;

        for (i = FLINT_MAX(band, 7); i <= rc; i++)
        {
            nn_ptr acc = LOGS_L(i);
            slong b = wp - i;

            flint_mpn_zero(acc, nc);
            if (2 * i + 1 <= wp)
            {
                /* ones on bits [64 nc - 2i - 1, 64 nc - i) */
                slong lo = wp - 2 * i - 1;
                slong j;
                for (j = lo; j < wp - i; j++)
                    acc[j / FLINT_BITS] |= UWORD(1) << (j % FLINT_BITS);
            }
            else
            {
                flint_mpn_store(acc, b / FLINT_BITS, ~UWORD(0));
                if (b % FLINT_BITS)
                    acc[b / FLINT_BITS] =
                        (UWORD(1) << (b % FLINT_BITS)) - 1;
            }
        }

        for (i = 7; i < FLINT_MIN(band, rc + 1); i++)
        {
            ulong g = _fixed_exp_logs[i * nc];
            if (g + FIXED_TAB_GUARD_SLACK < g)
                _fixed_tab_entry_exact(LOGS_L(i), 0, (ulong) i, nc);
        }
    }

}


/* Read-only view of the table for code outside the library (the
   thread-local storage itself is not exported; see fixed.h):
   the top n limbs of entry i, valid until the next ensure call on
   this thread. */
nn_srcptr
_fixed_exp_logs_entry(slong i, slong n)
{
    FLINT_ASSERT(i >= 0 && i <= _fixed_exp_logs_r);
    FLINT_ASSERT(n >= 1 && n <= _fixed_exp_logs_n - 1);
    return _fixed_exp_logs + i * _fixed_exp_logs_n
        + (_fixed_exp_logs_n - n);
}

/* free the cached table (tuning programs measure cold rebuilds;
   exported so that code outside the library never touches the
   thread-local data directly) */
void
_fixed_exp_logs_clear(void)
{
    flint_free(_fixed_exp_logs);
    _fixed_exp_logs = NULL;
    _fixed_exp_logs_n = 0;
    _fixed_exp_logs_r = 0;
}

/* the largest index the cached table currently covers */
slong
_fixed_exp_logs_max_index(void)
{
    return _fixed_exp_logs_r;
}




/* Product-form reconstruction for large ylen: accumulate
   P = prod (1 + 2^-used[j]) exactly in a growing lsb-anchored buffer
   (P = m 2^-B with B the accumulated bottom extent, so several
   factors fit in one limb, then two, three, ...), capping the extent
   at 64 ylen + 64 bits -- each capped factor then drops less than
   one bit below the cap, u such drops staying far under an output
   ulp -- and finish with a single multiplication y <- floor(y P)
   instead of one shift-and-add per factor.  With the used indices
   increasing, B grows only quadratically (about (sum of used
   indices) bits, r^2/4 in the worst case), so for large ylen the
   final multiplication is short times long and the whole
   reconstruction costs a fraction of the shift-and-add chain; the
   result is also tighter, one truncation instead of one per
   factor.  With the middle-product finish the path wins on the
   whole exp call from n ~ 256 at the default ladder on the
   development machine (7% there, 12-13% at n = 4096..8192, a wash
   only around n = 512 where the capped extent is at its relative
   largest).  Retune with the Newton cutoffs. */
static void
_fixed_exp_recon_prod(nn_ptr y, slong ylen, const slong * used,
    slong j, slong num)
{
    slong capd = FLINT_BITS * (ylen + 1);   /* keep bits above 2^-capd */
    slong B = 0, mlen = 1, alloc, plen;
    nn_ptr m, t, prod;
    TMP_INIT;

    TMP_START;
    /* the mantissa extent stays below capd + 63 + max_index bits:
       the depth cap plus one partial limb plus one factor's shift.
       A fixed "+ 16" slack here encoded max_index / 64 <= 12 and
       OVERFLOWED on 32-bit limbs, where the same indices span twice
       the limbs (an all-ones input at n = 420, r >= 576 corrupted
       the product a few limbs deep) */
    alloc = ylen + 4 + used[num - 1] / FLINT_BITS;
    m = TMP_ALLOC((2 * alloc + alloc + ylen + 2) * sizeof(ulong));
    t = m + alloc;
    prod = t + alloc;

    m[0] = 1;

    for (; j < num; j++)
    {
        slong i = used[j], q = i / FLINT_BITS, b = i % FLINT_BITS;
        slong tl, d;

        /* m' 2^-(B+i) = (m 2^i + m) 2^-(B+i) = m 2^-B (1 + 2^-i) */
        tl = mlen + q + (b != 0);
        flint_mpn_zero(t, q);
        if (b)
            t[mlen + q] = mpn_lshift(t + q, m, mlen, (int) b);
        else
            flint_mpn_copyi(t + q, m, mlen);
        while (tl > 1 && t[tl - 1] == 0)
            tl--;
        /* t = m 2^i (tl limbs, tl >= mlen); m' = t + m */
        {
            ulong cy = mpn_add(t, t, tl, m, mlen);
            t[tl] = cy;
            mlen = tl + (cy != 0);
        }
        { nn_ptr u_ = m; m = t; t = u_; }
        B += i;

        /* drop whole limbs sitting deeper than capd bits: each drop
           loses less than 2^-capd, one output guard bit, and there
           are at most r of them in total */
        d = (B - capd) / FLINT_BITS;
        if (d > 0)
        {
            flint_mpn_copyi(m, m + d, mlen - d);
            mlen -= d;
            B -= FLINT_BITS * d;
        }
    }

    /* y <- floor(y P) = (y m) >> B.  Only the product window at and
       above limb B/64 matters, so when the operands are nearly
       balanced a middle product beats the full multiplication.  Its
       lower-approximation deficit is bounded by
       min(ylen, mlen, zlo) 2^64 window units -- up to ~min ulps of
       the limb ABOVE the window bottom, a bound that all-ones
       operands very nearly attain -- so ONE sacrificial limb is not
       enough (at bit offset b = 0 up to min ulps would leak into
       the kept result); TWO are, leaving the deficit below
       min / 2^64 < 1 ulp of the result, hence within 1 ulp of the
       exact floor like the full product's plain truncation (the
       capped drops sit further below still), far inside the budget
       the per-factor chain used to spend. */
    plen = ylen + mlen;
    {
        slong q = B / FLINT_BITS, b = B % FLINT_BITS;

        if (2 * mlen >= ylen && q >= 2)
        {
            /* window two limbs below the kept result: both
               sacrificial limbs absorb the deficit */
            flint_mpn_mulmid(prod + q - 2, y, ylen, m, mlen,
                q - 2, plen);
            if (b)
                mpn_rshift(prod + q - 2, prod + q - 2,
                    plen - q + 2, (int) b);
            flint_mpn_copyi(y, prod + q, ylen);
        }
        else
        {
            if (mlen >= ylen)
                flint_mpn_mul(prod, m, mlen, y, ylen);
            else
                flint_mpn_mul(prod, y, ylen, m, mlen);
            if (b)
                mpn_rshift(prod + q, prod + q, plen - q, (int) b);
            flint_mpn_copyi(y, prod + q, ylen);
        }
    }
    TMP_END;
}

#ifndef FIXED_EXP_RECON_PROD_CUTOFF
#define FIXED_EXP_RECON_PROD_CUTOFF 192
#endif

/* multiply y (ylen limbs, fraction plus a units limb) by the factors
   (1 + 2^-used[j]) for j0 <= j < num: each is a shift and an add.
   The used indices are increasing; process them in chunks sharing the
   same limb offset.  The value stays below e, so no carry ever leaves
   the units limb.  sh is scratch of ylen limbs. */
void
_fixed_exp_recon(nn_ptr y, nn_ptr sh, slong ylen, const slong * used,
    slong j, slong num)
{
    if (ylen >= FIXED_EXP_RECON_PROD_CUTOFF && num - j >= 8)
    {
        _fixed_exp_recon_prod(y, ylen, used, j, num);
        return;
    }

    if (j == 0 && num > 0 && used[0] == 0)
    {
        mpn_add_n(y, y, y, ylen);               /* factor 2 */
        j = 1;
    }

    for (; j < num && used[j] < FLINT_BITS; j++)
    {
        mpn_rshift(sh, y, ylen, used[j]);
        mpn_add_n(y, y, sh, ylen);
    }

    while (j < num)
    {
        slong q = used[j] / FLINT_BITS;
        slong stop = (q + 1) * FLINT_BITS;
        slong len = ylen - q;

        for (; j < num && used[j] < stop; j++)
        {
            slong b = used[j] - q * FLINT_BITS;
            ulong cy;

            if (b != 0)
                mpn_rshift(sh, y + q, len, b);
            else
                flint_mpn_copyi(sh, y + q, len);

            cy = mpn_add_n(y, y, sh, len);
            mpn_add_1(y + len, y + len, q, cy);
        }
    }
}


/* Greedy bitwise reduction of t (wn limbs) by L_0, ..., L_r,
   recording the used indices; returns their number.  The subtractions
   are decided one FLINT_BITS window at a time on a single-limb model
   h of the remainder: h is kept a lower bound on the remainder's
   current window limb, with upper slack e growing by one per model
   subtraction (each table limb is a truncation), so h > lt proves
   L_i <= t and h <= lt - e - 1 proves t < L_i.  Decisions inside the
   ambiguity band [lt - e, lt] -- a vanishing fraction of random
   steps, plus the last few steps of each window where lt itself is
   small -- fall back to the exact compare-subtract, as does the first
   step of each window (the remainder may still carry one bit above
   it).  Decided subtractions are recorded and applied to the
   full-precision remainder in unconditional batches, so the bulk of
   the loop touches one limb per step and the full-width work is
   branch-free.  The decisions agree exactly with the plain
   full-precision greedy. */
/* Shared greedy reduction: subtract tab[i] (entries of tabn limbs,
   value tab_i < 2^-i with tab_{i-1} < 2 tab_i) from (t, wn) for
   i = istart..r whenever it fits, recording the used indices; used
   by exp (tab_i = log(1 + 2^-i), istart = 0) and by the rotation
   reduction of sin/cos (tab_i = atan(2^-i), istart = 1, applied to
   half the angle -- the entries MUST satisfy tab_i < 2^-i, which
   the doubled angles 2 atan(2^-i) ~ 2^(1-i) would not). */
slong
_fixed_bitwise_reduce(nn_ptr t, slong wn, int r, slong istart,
    nn_srcptr tab, slong tabn, slong * used)
{
    slong num = 0, bj = 0, nc = tabn;
    slong i, c;

    for (c = 0; FLINT_BITS * c <= r; c++)
    {
        slong i0 = FLINT_MAX(FLINT_BITS * c, istart);
        slong i1 = FLINT_MIN((slong) r, FLINT_BITS * (c + 1) - 1);
        ulong h, e;
        nn_srcptr Lp;

        /* apply pending subtractions from the previous window */
        for (; bj < num; bj++)
            mpn_sub_n(t, t, TAB_E(used[bj]) + (nc - wn), wn);

        /* Exact step at the window boundary.  This must leave
           t < tab[i0] < 2^-i0, which is precisely the precondition
           of the single-limb model below (that t[wn - 1 - c] is the
           leading limb).  A single subtraction is not enough to
           guarantee it: the table entries are floors, so each
           accepted index under-subtracts by up to one ulp, and this
           creep can push t above tab[i0 - 1] -- close to 2 tab[i0]
           -- whenever the entry's deficit from 2^-i0 is itself
           below an ulp.  Repeating the subtraction restores the
           invariant unconditionally; an index used twice simply
           means its factor is applied twice, which the callers'
           reconstructions handle (the decomposition stays bounded by
           the input, not by the table). */
        Lp = TAB_E(i0) + (nc - wn);
        while (mpn_cmp(t, Lp, wn) >= 0)
        {
            mpn_sub_n(t, t, Lp, wn);
            used[num++] = i0;
        }
        bj = num;

        h = t[wn - 1 - c];
        e = 0;

        for (i = i0 + 1; i <= i1; i++)
        {
            ulong lt = tab[i * nc + (nc - 1 - c)];

            if (h - lt + e <= 2 * e)
            {
                /* ambiguity band: apply pending and decide exactly */
                for (; bj < num; bj++)
                    mpn_sub_n(t, t, TAB_E(used[bj]) + (nc - wn), wn);
                Lp = TAB_E(i) + (nc - wn);
                if (mpn_cmp(t, Lp, wn) >= 0)
                {
                    mpn_sub_n(t, t, Lp, wn);
                    used[num++] = i;
                }
                bj = num;
                h = t[wn - 1 - c];
                e = 0;
            }
            else
            {
                ulong sub = (h > lt);

                h -= (lt + 1) & (-sub);
                used[num] = i;
                num += sub;
                e += sub;
            }
        }
    }

    for (; bj < num; bj++)
        mpn_sub_n(t, t, TAB_E(used[bj]) + (nc - wn), wn);

    /* the truncated table lets the remainder creep marginally above
       tab[r]; extra subtractions restore t < tab[r] < 2^-r */
    {
        nn_srcptr Lr = TAB_E(r) + (nc - wn);

        while (mpn_cmp(t, Lr, wn) >= 0)
        {
            mpn_sub_n(t, t, Lr, wn);
            used[num++] = r;
        }
    }

    return num;
}

#if FLINT_BITS == 64
#endif

/* smallest n at which each reduction parameter becomes optimal
   (generated by src/fixed/tune/tune-bitwise-r.c with nmax deep enough
   to reach the r = 768 end of the ladder; x86-64 defaults) */
static const int _fixed_exp_bitwise_rs_r_tab[] =
    {32, 64, 128, 192, 256, 320, 384, 448, 512, 576, 640, 704, 768};
static const short _fixed_exp_bitwise_rs_n_tab[] =
    {1, 12, 41, 71, 100, 182, 258, 282, 337, 414, 486, 496, 567};

int
fixed_exp_bitwise_rs_default_r(slong n)
{
    slong j;

#if FLINT_BITS == 64
    /* the compile-time constants of exp_opt_<n>.c: keep in sync with
       the dev/tune_fixed.py --pin values that emitted those files */
    static const int opt_r[] = {0, 12, 16, 16, 16, 16, 24, 32};
    if (n <= 7)
        return opt_r[n];
#endif

    for (j = 0; j + 1 < (slong) (sizeof(_fixed_exp_bitwise_rs_n_tab)
            / sizeof(short))
            && n >= _fixed_exp_bitwise_rs_n_tab[j + 1]; j++)
        ;
    return _fixed_exp_bitwise_rs_r_tab[j];
}

void
fixed_exp_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)
{
#if FLINT_BITS == 64
    int r0;
#endif
    slong num, wn;
    nn_ptr t, y, sh;
    slong * used;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);
    /* on 32-bit limbs, n = 1 would clamp r to 16 below, violating
       the 2^-32 contract of the series functions */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

#if FLINT_BITS == 64
    r0 = r;         /* the specialized sizes trigger on r = 0 only */
#endif

    if (r == 0)
        r = fixed_exp_bitwise_rs_default_r(n);

#if FLINT_BITS == 64
    /* the fully specialized per-size implementations (exp_opt_<n>.c,
       emitted and tuned by dev/tune_fixed.py) serve the default
       dispatch: with r = 0 the caller has left the reduction
       parameter to us, and the whole call collapses to straight-line
       code with a compile-time r and the series built for it */
    if (r0 == 0 && n <= 7)
    {
        static void (* const tab[])(nn_ptr, nn_srcptr) = {
            NULL, fixed_exp_opt_1, fixed_exp_opt_2, fixed_exp_opt_3,
            fixed_exp_opt_4, fixed_exp_opt_5, fixed_exp_opt_6,
            fixed_exp_opt_7
        };
        tab[n](res, x);
        return;
    }
#endif

    r = FLINT_MAX(r, 32);
    r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);

    wn = n;

    _fixed_exp_logs_ensure(wn, r);

    TMP_START;
    t = TMP_ALLOC((wn + 2 * (n + 1)) * sizeof(ulong));
    y = t + wn;
    sh = y + (n + 1);
    used = TMP_ALLOC(FIXED_BITWISE_REDUCE_USED_ALLOC(r) * sizeof(slong));

    flint_mpn_copyi(t, x, n);

    num = _fixed_bitwise_reduce(t, wn, r, 0, _fixed_exp_logs,
        _fixed_exp_logs_n, used);

    /* exp of the reduced argument */
    fixed_exp_reduced(y, t, wn, (flint_bitcnt_t) r, 0);

    _fixed_exp_recon(y, sh, n + 1, used, 0, num);

    flint_mpn_copyi(res, y, n + 1);

    TMP_END;
}
