/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "padic_radix.h"
#include "arb/impl.h"
#include "gr.h"

/* For 2 <= b - a < this bound, and when the block value x = p^vr ru fits in a
   single limb, the binary-splitting leaf is evaluated by an iterative mpn
   accumulation rather than by recursion.  Tuning parameter. */
#define PADIC_RADIX_LOG_BSPLIT_BASECASE 64
#define PADIC_RADIX_LOG_BC_LIMBS (2 * PADIC_RADIX_LOG_BSPLIT_BASECASE + 4)

/* Digit width of the first (lowest-valuation) block peeled from the argument.
   Tuning parameter; must be a power of two. */
#define PADIC_RADIX_LOG_SMALLEST_BLOCK 4

/*
    Iterative basecase for a block whose value x = p^vr ru fits in a single limb
    xl = p^vr ru.  Accumulates the binary-splitting leaf (Tu, B) for [a, b) in
    plain mpn limbs and converts to radix once at the end.

    For the series sum_{i=a}^{b-1} x^i / i (local powers), T/B = sum with B the
    product a(a+1)...(b-1).  Tracking BXp = B * x^m (m = number of terms so far)
    keeps every multiplier single-limb:

        T   <- J*T + BXp*x,   B <- J*B,   BXp <- BXp*J*x   (J = a+1, ..., b-1),

    started from T = x, B = a, BXp = a*x.  The accumulated T is the full
    numerator p^vr Tu, so Tu is recovered by a final shift down by vr digits.
*/
static void
_padic_radix_log_bsplit_basecase_mpn(radix_integer_t Tu, radix_integer_t B,
    slong a, slong b, ulong xl, slong vr, slong l, const radix_t radix)
{
    ulong T[PADIC_RADIX_LOG_BC_LIMBS];
    ulong BB[PADIC_RADIX_LOG_BC_LIMBS];
    ulong BXp[PADIC_RADIX_LOG_BC_LIMBS];
    slong Tn, Bn, BXpn, J, need;
    ulong hi, lo, cy;
    nn_ptr d;

    /* index a:  T = x,  B = a,  BXp = a*x */
    T[0] = xl;          Tn = 1;
    BB[0] = (ulong) a;  Bn = 1;
    umul_ppmm(hi, lo, (ulong) a, xl);
    BXp[0] = lo;
    BXp[1] = hi;
    BXpn = 1 + (hi != 0);

    for (J = a + 1; J < b; J++)
    {
        /* T = J*T + BXp*x:  scale T by J, then fuse the BXp*x add in one pass
           over the (usually equal-length) overlap with mpn_addmul_1. */
        hi = mpn_mul_1(T, T, Tn, (ulong) J);
        T[Tn] = hi;  Tn += (hi != 0);

        if (Tn < BXpn)
        {
            flint_mpn_zero(T + Tn, BXpn - Tn);
            Tn = BXpn;
        }
        cy = mpn_addmul_1(T, BXp, BXpn, xl);
        if (cy != 0)
        {
            if (Tn > BXpn)
                cy = mpn_add_1(T + BXpn, T + BXpn, Tn - BXpn, cy);
            if (cy != 0)
            {
                T[Tn] = cy;  Tn++;
            }
        }

        /* B = J*B */
        hi = mpn_mul_1(BB, BB, Bn, (ulong) J);
        BB[Bn] = hi;  Bn += (hi != 0);

        /* BXp = BXp*J*x  (skip the final, unused update); one mul_1 if J*x fits */
        if (J + 1 < b)
        {
            umul_ppmm(hi, lo, (ulong) J, xl);
            if (hi == 0)
            {
                hi = mpn_mul_1(BXp, BXp, BXpn, lo);
                BXp[BXpn] = hi;  BXpn += (hi != 0);
            }
            else
            {
                hi = mpn_mul_1(BXp, BXp, BXpn, (ulong) J);
                BXp[BXpn] = hi;  BXpn += (hi != 0);
                hi = mpn_mul_1(BXp, BXp, BXpn, xl);
                BXp[BXpn] = hi;  BXpn += (hi != 0);
            }
        }
    }

    /* Tu = T / p^vr (T has valuation vr in p-adic digits, vr < e). */
    need = radix_set_mpn_need_alloc(Tn, radix);
    d = radix_integer_fit_limbs(Tu, need, radix);
    Tu->size = radix_set_mpn(d, T, Tn, radix);
    if (vr > 0)
        radix_integer_rshift_digits(Tu, Tu, vr, radix);
    radix_integer_mod_limbs(Tu, Tu, l, radix);

    need = radix_set_mpn_need_alloc(Bn, radix);
    d = radix_integer_fit_limbs(B, need, radix);
    B->size = radix_set_mpn(d, BB, Bn, radix);
    radix_integer_mod_limbs(B, B, l, radix);
}

/*
    Binary splitting of  sum_{i=a}^{b-1} x^i / i  for a block x = p^vr * ru,
    where ru is a unit.  The factor of p^vr carried by every power of x is kept
    out of the multiplications: the series numerator T = p^vr * Tu is returned
    through its unit-scaled part Tu, and B = (b-1)!/(a-1)! is the denominator.
    Powers ru^step are read from the precomputed table xpow (indexed through the
    bs-exponent table xexp); the p^vr factors enter as digit shifts.  All values
    are modulo B^l.  Identities, with T = p^vr Tu and S(a,b) = T/B:

        Tu(a,a+1) = ru,                         B(a,a+1) = a
        Tu(a,a+2) = (a+1) ru + a p^vr ru^2,     B(a,a+2) = a(a+1)
        Tu        = Tu_L B_R + p^{vr step} ru^step B_L Tu_R,   step = (b-a)/2
        B         = B_L B_R
*/
static void
_padic_radix_log_bsplit(radix_integer_t Tu, radix_integer_t B,
    slong a, slong b, const slong * xexp, const radix_integer_struct * xpow,
    slong vr, slong l, ulong xl, int xl_fits, const radix_t radix)
{
    slong e = radix->exp;
    slong el = e * l;

    if (xl_fits && (b - a) < PADIC_RADIX_LOG_BSPLIT_BASECASE)
    {
        /* 1 <= b - a < tuning, block fits a limb: iterative mpn leaf.  This
           also covers b - a == 1 (empty loop, Tu = ru), so the power table is
           never consulted along this path. */
        _padic_radix_log_bsplit_basecase_mpn(Tu, B, a, b, xl, vr, l, radix);
    }
    else if (b - a == 1)
    {
        radix_integer_set(Tu, xpow + 0, radix);          /* ru^1 */
        radix_integer_set_ui(B, (ulong) a, radix);
    }
    else if (b - a == 2)
    {
        slong i2 = _arb_get_exp_pos(xexp, 2);
        radix_integer_t tmp;
        radix_integer_init(tmp, radix);

        /* Tu = (a+1) ru + a p^vr ru^2 */
        radix_integer_set_ui(tmp, (ulong) (a + 1), radix);
        radix_integer_mullow_limbs(Tu, xpow + 0, tmp, l, radix);
        if (vr < el)
        {
            radix_integer_set_ui(tmp, (ulong) a, radix);
            radix_integer_mullow_limbs(tmp, xpow + i2, tmp, l, radix);  /* a ru^2 */
            radix_integer_lshift_digits(tmp, tmp, vr, radix);
            radix_integer_mod_limbs(tmp, tmp, l, radix);
            radix_integer_add(Tu, Tu, tmp, radix);
            radix_integer_mod_limbs(Tu, Tu, l, radix);
        }

        /* B = a(a+1) */
        radix_integer_set_ui(B, (ulong) a, radix);
        radix_integer_set_ui(tmp, (ulong) (a + 1), radix);
        radix_integer_mullow_limbs(B, B, tmp, l, radix);

        radix_integer_clear(tmp, radix);
    }
    else
    {
        slong step = (b - a) / 2;
        slong m = a + step;
        slong shift;
        int include;
        radix_integer_t TuR, BR, t2;

        radix_integer_init(TuR, radix);
        radix_integer_init(BR, radix);
        radix_integer_init(t2, radix);

        _padic_radix_log_bsplit(Tu, B, a, m, xexp, xpow, vr, l, xl, xl_fits, radix);
        _padic_radix_log_bsplit(TuR, BR, m, b, xexp, xpow, vr, l, xl, xl_fits, radix);

        /* Tu = Tu_L B_R + p^{vr step} ru^step B_L Tu_R */
        radix_integer_mullow_limbs(Tu, Tu, BR, l, radix);   /* Tu_L B_R */

        if (step != 0 && vr > el / step)
            include = 0;
        else
        {
            shift = vr * step;
            include = (shift < el);
        }

        if (include)
        {
            slong is = _arb_get_exp_pos(xexp, step);
            radix_integer_mullow_limbs(t2, xpow + is, B, l, radix);   /* ru^step B_L */
            radix_integer_mullow_limbs(t2, t2, TuR, l, radix);        /* * Tu_R */
            radix_integer_lshift_digits(t2, t2, shift, radix);
            radix_integer_mod_limbs(t2, t2, l, radix);
            radix_integer_add(Tu, Tu, t2, radix);
            radix_integer_mod_limbs(Tu, Tu, l, radix);
        }

        radix_integer_mullow_limbs(B, B, BR, l, radix);     /* B_L B_R */

        radix_integer_clear(TuR, radix);
        radix_integer_clear(BR, radix);
        radix_integer_clear(t2, radix);
    }
}

/*
    rop = sum_{i>=1} x^i / i  modulo p^N, the negative of log(1 - x), for a
    block x of p-adic valuation at least w (1 <= w < N).  The result is the
    nonnegative residue and carries the valuation vr = v_p(x): only the unit
    part ru = x / p^vr enters the multiplications, and the result is shifted up
    by vr at the end, so the series is evaluated to the reduced precision
    N - vr.
*/
static void
_padic_radix_log_bsplit_block(radix_integer_t rop, const radix_integer_t x,
    slong w, slong N, const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    slong e = radix->exp;
    slong vr, Nr, n, k, l, lNr, kk, length, i;
    slong * xexp;
    radix_integer_struct * xpow;
    radix_integer_t ru, Tu, B, ub;
    ulong xl = 0, ru_l;
    int xl_fits, use_table;

    vr = radix_integer_valuation_digits(x, radix);
    if (vr >= N)                                  /* contribution below p^N */
    {
        radix_integer_zero(rop, radix);
        return;
    }

    Nr = N - vr;                                  /* relevant precision */

    n = _padic_radix_log_bound(FLINT_MAX(w, vr), N, p);
    n = FLINT_MAX(n, 2);

    k = (n >= 2) ? (n - 2) / (slong) (p - 1) : 0; /* guard >= v_p((n-1)!) */
    lNr = (Nr + e - 1) / e;
    l = lNr + (k + e - 1) / e;                    /* room for the p^kk strip */

    radix_integer_init(ru, radix);
    radix_integer_init(Tu, radix);
    radix_integer_init(B, radix);
    radix_integer_init(ub, radix);

    radix_integer_rshift_digits(ru, x, vr, radix);    /* unit part */
    radix_integer_mod_limbs(ru, ru, l, radix);

    /* The block value p^vr ru fits in a single limb exactly when vr < e and the
       unit is below p^{e-vr}; the iterative mpn leaf then applies. */
    ru_l = (ru->size == 0) ? 0 : ru->d[0];
    xl_fits = (FLINT_ABS(ru->size) <= 1) && (vr < e) && (ru_l < radix->bpow[e - vr]);
    if (xl_fits)
        xl = radix->bpow[vr] * ru_l;              /* block value, < B */

    /* The power table is only consulted by the recursive/leaf paths; if the
       whole series [1, n) is taken by the iterative basecase, skip building it. */
    use_table = !(xl_fits && (n - 1) < PADIC_RADIX_LOG_BSPLIT_BASECASE);

    xexp = NULL;
    xpow = NULL;
    length = 0;

    if (use_table)
    {
        /* Power table ru^xexp[i] for the binary splitting over [1, n). */
        xexp = flint_calloc(2 * FLINT_BITS, sizeof(slong));
        length = _arb_compute_bs_exponents(xexp, n - 1);
        xpow = flint_malloc(sizeof(radix_integer_struct) * length);
        for (i = 0; i < length; i++)
            radix_integer_init(xpow + i, radix);

        radix_integer_set(xpow + 0, ru, radix);           /* ru^1 */
        for (i = 1; i < length; i++)
        {
            if (xexp[i] == 2 * xexp[i - 1])
                radix_integer_mullow_limbs(xpow + i, xpow + i - 1, xpow + i - 1, l, radix);
            else if (xexp[i] == 2 * xexp[i - 2])
                radix_integer_mullow_limbs(xpow + i, xpow + i - 2, xpow + i - 2, l, radix);
            else if (xexp[i] == 2 * xexp[i - 1] + 1)
            {
                radix_integer_mullow_limbs(xpow + i, xpow + i - 1, xpow + i - 1, l, radix);
                radix_integer_mullow_limbs(xpow + i, xpow + i, ru, l, radix);
            }
            else if (xexp[i] == 2 * xexp[i - 2] + 1)
            {
                radix_integer_mullow_limbs(xpow + i, xpow + i - 2, xpow + i - 2, l, radix);
                radix_integer_mullow_limbs(xpow + i, xpow + i, ru, l, radix);
            }
            else
                flint_throw(FLINT_ERROR, "padic_radix log: power table malformed\n");
        }
    }

    _padic_radix_log_bsplit(Tu, B, 1, n, xexp, xpow, vr, l, xl, xl_fits, radix);

    /* Strip the common power of p (= v_p(B) = v_p(Tu)) so both are units. */
    kk = radix_integer_valuation_digits(B, radix);
    if (kk > 0)
    {
        radix_integer_rshift_digits(B, B, kk, radix);
        radix_integer_rshift_digits(Tu, Tu, kk, radix);
    }

    /* ub = Tu / B  (mod p^Nr); rop = p^vr * ub. */
    if (!radix_integer_divmod_limbs(ub, Tu, B, lNr, radix))
        flint_throw(FLINT_ERROR, "_padic_radix_log_bsplit_block: series "
            "denominator is not invertible after stripping\n");

    radix_integer_mod_digits(ub, ub, Nr, radix);
    radix_integer_lshift_digits(rop, ub, vr, radix);
    radix_integer_mod_digits(rop, rop, N, radix);

    if (use_table)
    {
        for (i = 0; i < length; i++)
            radix_integer_clear(xpow + i, radix);
        flint_free(xpow);
        flint_free(xexp);
    }
    radix_integer_clear(ru, radix);
    radix_integer_clear(Tu, radix);
    radix_integer_clear(B, radix);
    radix_integer_clear(ub, radix);
}

/*
    rop = log(1 - y) modulo p^N, with y treated as an exact integer of p-adic
    valuation at least the convergence threshold and less than N.

    Uses the balanced decomposition 1 - y = prod_j (1 - r_j), where each r_j is
    a doubling-width block peeled off the running cofactor t and the cofactor is
    corrected by (1 - r_j)^{-1}; then log(1 - y) = sum_j log(1 - r_j).  The
    nonnegative residue is returned.  Always returns GR_SUCCESS (an internal
    invariant failure aborts via flint_throw).
*/
void
_padic_radix_log_balanced(radix_integer_t rop, const radix_integer_t y, slong N,
    const radix_t radix)
{
    slong e = radix->exp;
    slong lN, bound, w;
    radix_integer_t t, r, om, blk, z;

    if (N <= 0)
    {
        radix_integer_zero(rop, radix);
        return;
    }

    lN = (N + e - 1) / e;

    radix_integer_init(t, radix);
    radix_integer_init(r, radix);
    radix_integer_init(om, radix);
    radix_integer_init(blk, radix);
    radix_integer_init(z, radix);

    radix_integer_mod_digits(t, y, N, radix);     /* t = y mod p^N */
    radix_integer_zero(z, radix);

    bound = PADIC_RADIX_LOG_SMALLEST_BLOCK;   /* truncation p^bound of first block */
    w = 1;                                     /* valuation lower bound of first block */
    while (!radix_integer_is_zero(t, radix))
    {
        radix_integer_mod_digits(r, t, bound, radix);   /* low block */
        radix_integer_sub(t, t, r, radix);              /* zero low `bound` digits */

        if (!radix_integer_is_zero(t, radix))
        {
            /* t <- (t - r) / (1 - r)  (mod p^N); 1 - r is a 1-unit */
            radix_integer_one(om, radix);
            radix_integer_sub(om, om, r, radix);
            if (!radix_integer_divmod_limbs(t, t, om, lN, radix))
                flint_throw(FLINT_ERROR, "_padic_radix_log_balanced: 1 - r is "
                    "not invertible\n");
            radix_integer_mod_digits(t, t, N, radix);
        }

        if (!radix_integer_is_zero(r, radix))
        {
            _padic_radix_log_bsplit_block(blk, r, w, N, radix);   /* sum r^i/i */
            radix_integer_sub(z, z, blk, radix);            /* += log(1 - r) */
        }

        w = bound;
        bound *= 2;
    }

    radix_integer_mod_digits(rop, z, N, radix);

    radix_integer_clear(t, radix);
    radix_integer_clear(r, radix);
    radix_integer_clear(om, radix);
    radix_integer_clear(blk, radix);
    radix_integer_clear(z, radix);
}
