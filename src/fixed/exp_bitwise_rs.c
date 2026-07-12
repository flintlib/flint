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
#include "impl.h"

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

#include "frac_bsplit.inc"

#define LOGS_L(i) (_fixed_exp_logs + (i) * _fixed_exp_logs_n)
#define TAB_E(i) (tab + (i) * nc)

/* floor(x 2^(FLINT_BITS nc)) of a nonnegative arb into nc limbs */
static void
_fixed_tab_store(nn_ptr e, const arb_t x, slong nc)
{
    arf_t lb;
    fmpz_t f;

    arf_init(lb);
    fmpz_init(f);
    arb_get_lbound_arf(lb, x, FLINT_BITS * nc + 30);
    arf_mul_2exp_si(lb, lb, FLINT_BITS * nc);
    arf_get_fmpz(f, lb, ARF_RND_FLOOR);
    FLINT_ASSERT(fmpz_sgn(f) >= 0);
    fmpz_get_ui_array(e, nc, f);
    arf_clear(lb);
    fmpz_clear(f);
}

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
   all of them; beyond that, log(1 + 2^-i) = 2 atanh(1/(2^(i+1) + 1))
   through arb_atan1_frac_bsplit (frac_bsplit.inc), switching to the
   direct log(1 + 2^-i) series of arb_log1_frac_bsplit once the
   precision passes 6000 bits and i passes 30, where its plainer
   recursion wins.  TODO: reimplement the binary splitting natively
   with mpn arithmetic; arb is fine for now, and this tier is not
   the bottleneck.

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
                    if (prec >= 6000 && i >= 30)
                    {
                        /* direct log(1 + 1/2^i) series */
                        fmpz_one(p);
                        fmpz_one(q);
                        fmpz_mul_2exp(q, q, i);
                        arb_log1_frac_bsplit(x, p, q, wp);
                    }
                    else
                    {
                        /* log(1 + 2^-i) = 2 atanh(1 / (2^(i+1) + 1)) */
                        fmpz_one(p);
                        fmpz_one(q);
                        fmpz_mul_2exp(q, q, i + 1);
                        fmpz_add_ui(q, q, 1);
                        arb_atan1_frac_bsplit(x, p, q, 1, wp);
                        arb_mul_2exp_si(x, x, 1);
                    }
                    break;
            }

            _fixed_tab_store(LOGS_L(i), x, nc);
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
}


/* Read-only view of the table for code outside the library (the
   thread-local storage itself is not exported; see impl.h):
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

/* the largest index the cached table currently covers */
slong
_fixed_exp_logs_max_index(void)
{
    return _fixed_exp_logs_r;
}


/* Use the sinh series (half the terms) plus a squaring and a square
   root once the direct exp series gets long enough.  Measured
   crossovers on x86-64: for r < 64 (series through the pre32 path)
   sinh wins from about 45 terms (n ~ 24 at r = 32); for r >= 64 the
   windowed pre64 series is more efficient per term and mpn_sqrtrem
   sets a higher floor, moving the crossover to about 128 terms
   (n ~ 2r).  The margins near the crossovers are within a few
   percent, so the exact placement is not critical. */
#define EXP_USE_SINH(wn, r) \
    (FLINT_BITS * (wn) >= (((r) >= 64) ? 128 : 45) * (slong) (r))

/* exp(t) of the reduced argument t < 2^-r into y (wn + 1 limbs:
   wn fraction limbs and a units limb); the series functions pick the
   32- or 64-bit internal range from the top limb of t.  When the
   direct series would need many terms, evaluate sinh(t) instead --
   the odd series has half the terms -- and reconstruct
   exp(t) = sinh(t) + sqrt(1 + sinh(t)^2), costing one squaring and
   one square root.  The sinh error (FIXED_SINH_RS_MAX_ERR = 15 ulp),
   the truncated squaring and the floored integer square root
   together contribute some 20 ulp, amplified below e by the
   reconstruction: this is part of the constant term of
   FIXED_EXP_BITWISE_RS_MAX_ERR. */
static void
_fixed_exp_reduced(nn_ptr y, nn_srcptr t, slong wn, int r, int use_sinh)
{
    if (!use_sinh)
    {
        fixed_exp_rs(y, t, wn);
    }
    else
    {
        nn_ptr s, u2, rt, rem;
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(((wn + 1) + (2 * wn + 1) + (wn + 1) + (wn + 2))
            * sizeof(ulong));
        u2 = s + (wn + 1);
        rt = u2 + (2 * wn + 1);
        rem = rt + (wn + 1);

        fixed_sinh_rs(s, t, wn);

        /* u2 = (1 + sinh(t)^2) 2^(2 FLINT_BITS wn) */
        flint_mpn_zero(u2, wn);
        flint_mpn_sqrhigh(u2 + wn, s, wn);
        u2[2 * wn] = 1;

        /* cosh(t) = sqrt(1 + sinh(t)^2) */
        mpn_sqrtrem(rt, rem, u2, 2 * wn + 1);

        /* exp(t) = sinh(t) + cosh(t) */
        mpn_add_n(y, rt, s, wn + 1);

        TMP_END;
    }
}


/* multiply y (ylen limbs, fraction plus a units limb) by the factors
   (1 + 2^-used[j]) for j0 <= j < num: each is a shift and an add.
   The used indices are increasing; process them in chunks sharing the
   same limb offset.  The value stays below e, so no carry ever leaves
   the units limb.  sh is scratch of ylen limbs. */
void
_fixed_exp_recon(nn_ptr y, nn_ptr sh, slong ylen, const slong * used,
    slong j, slong num)
{
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
    _fixed_exp_reduced(y, t, wn, r, EXP_USE_SINH(wn, r));

    _fixed_exp_recon(y, sh, n + 1, used, 0, num);

    flint_mpn_copyi(res, y, n + 1);

    TMP_END;
}
