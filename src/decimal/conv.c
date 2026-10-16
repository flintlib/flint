/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "qqbar.h"
#include "decimal.h"
#include "double_extras.h"
#include "mag.h"
#include "fmpq.h"
#include "fmpz_extras.h"
#include "arf.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "mpn_extras.h"

/* ------------------------------------------------------------------------- */
/*    From integers                                                          */
/* ------------------------------------------------------------------------- */

/* Round the (signed) mantissa currently stored in res->m with limb
   exponent exp (fmpz, may be zero). */
static int
_decfloat_round_own_mantissa(decfloat_t res, const fmpz_t exp, slong prec, int rnd, gr_ctx_t ctx)
{
    slong n = FLINT_ABS(res->m.size);
    int negative = res->m.size < 0;

    if (n == 0)
        return decfloat_zero(res, ctx);

    radix_integer_fit_limbs(&res->m, n + 1, DECIMAL_CTX_RADIX(ctx));
    return _decfloat_set_round_limbs(res, res->m.d, n, negative, exp, 0, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_set_round_fmpz(decfloat_t res, const fmpz_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz zero = 0;

    if (fmpz_is_zero(x))
        return decfloat_zero(res, ctx);

    radix_integer_set_fmpz(&res->m, x, DECIMAL_CTX_RADIX(ctx));
    return _decfloat_round_own_mantissa(res, &zero, prec, rnd, ctx);
}

int
decfloat_set_round_ui(decfloat_t res, ulong x, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz zero = 0;

    if (x == 0)
        return decfloat_zero(res, ctx);

    radix_integer_set_ui(&res->m, x, DECIMAL_CTX_RADIX(ctx));
    return _decfloat_round_own_mantissa(res, &zero, prec, rnd, ctx);
}

int
decfloat_set_round_si(decfloat_t res, slong x, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz zero = 0;

    if (x == 0)
        return decfloat_zero(res, ctx);

    radix_integer_set_si(&res->m, x, DECIMAL_CTX_RADIX(ctx));
    return _decfloat_round_own_mantissa(res, &zero, prec, rnd, ctx);
}

int
decfloat_set_round_fmpz_10exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t t, slong prec, int rnd, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    fmpz_t q, r, mm;
    ulong rr;
    int status;

    if (fmpz_is_zero(m))
        return decfloat_zero(res, ctx);

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(mm);

    fmpz_fdiv_q_ui(q, t, e);        /* q = floor(t/e) */
    fmpz_mul_ui(r, q, e);
    fmpz_sub(r, t, r);              /* r = t - q e in [0, e) */
    rr = fmpz_get_ui(r);

    if (rr == 0)
        fmpz_set(mm, m);
    else
        fmpz_mul_ui(mm, m, DECIMAL_CTX_RADIX(ctx)->bpow[rr]);

    radix_integer_set_fmpz(&res->m, mm, DECIMAL_CTX_RADIX(ctx));
    status = _decfloat_round_own_mantissa(res, q, prec, rnd, ctx);

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(mm);

    return status;
}

int
decfloat_set_fmpz_10exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
{
    return decfloat_set_round_fmpz_10exp_fmpz(res, m, e, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_set_fmpz_10exp_si(decfloat_t res, const fmpz_t m, slong e, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, e);
    status = decfloat_set_fmpz_10exp_fmpz(res, m, t, ctx);
    fmpz_clear(t);
    return status;
}

int
decfloat_set_si_10exp_si(decfloat_t res, slong m, slong e, gr_ctx_t ctx)
{
    fmpz_t t, u;
    int status;
    fmpz_init_set_si(t, e);
    fmpz_init_set_si(u, m);
    status = decfloat_set_fmpz_10exp_fmpz(res, u, t, ctx);
    fmpz_clear(t);
    fmpz_clear(u);
    return status;
}

/* Size limit (in bits) for the integers formed in the exact conversion of
   m 2^t to prec digits; beyond this, scaling in arb is cheaper. */
#define DECIMAL_EXACT_CONV_LIMIT(prec) (16 * _decimal_digits_to_bits(prec) + 4096)

/* 5^k for 0 <= k <= DECIMAL_POW5_MAX_EXP */
static const ulong _pow5_tab[DECIMAL_POW5_MAX_EXP + 1] = {
    UWORD(1), UWORD(5), UWORD(25), UWORD(125), UWORD(625), UWORD(3125),
    UWORD(15625), UWORD(78125), UWORD(390625), UWORD(1953125), UWORD(9765625),
    UWORD(48828125), UWORD(244140625), UWORD(1220703125),
#if FLINT_BITS == 64
    UWORD(6103515625), UWORD(30517578125), UWORD(152587890625), UWORD(762939453125),
    UWORD(3814697265625), UWORD(19073486328125), UWORD(95367431640625),
    UWORD(476837158203125), UWORD(2384185791015625), UWORD(11920928955078125),
    UWORD(59604644775390625), UWORD(298023223876953125),
    UWORD(1490116119384765625), UWORD(7450580596923828125)
#endif
};

/* Fast path for the scaled conversion: |m| fits in SMALL_MN limbs, 0 <= s
   (at most 12 DECIMAL_POW5_MAX_EXP) and the scaled integer fits in
   SMALL_LIMBS limbs (a few hundred digits). Computes q = floor(|m| 5^s
   2^(t+s)) with a sticky flag and rounds, using stack buffers. */
#define SMALL_LIMBS 16
#define SMALL_MN 8

/* Sets q (with room for SMALL_LIMBS + 1 limbs) to floor(m 5^s 2^(t+s))
   where m = {mp, mn} with mn <= SMALL_MN and 0 <= s, and sticky to whether the
   floor is inexact. Returns the length of q (0 if q = 0), or -1 if the
   sizes are not small. */
static slong
_scaled_floor(nn_ptr q, int * sticky, nn_srcptr mp, slong mn, slong t, slong s)
{
    slong qn, sh, bits, i, ss;

    if (s < 0 || s > 12 * DECIMAL_POW5_MAX_EXP || mn > SMALL_MN)
        return -1;

    flint_mpn_copyi(q, mp, mn);
    qn = mn;

    /* q = |m| 5^s, in chunks of 5^DECIMAL_POW5_MAX_EXP */
    for (ss = s; ss > 0; )
    {
        slong c = FLINT_MIN(ss, DECIMAL_POW5_MAX_EXP);
        if (qn == SMALL_LIMBS)
            return -1;
        q[qn] = mpn_mul_1(q, q, qn, _pow5_tab[c]);
        qn += (q[qn] != 0);
        ss -= c;
    }

    sh = t + s;
    *sticky = 0;

    if (sh > 0)
    {
        slong limbs = sh / FLINT_BITS, r = sh % FLINT_BITS;

        bits = (qn - 1) * FLINT_BITS + FLINT_BIT_COUNT(q[qn - 1]);
        if (bits + sh > SMALL_LIMBS * FLINT_BITS)
            return -1;

        if (limbs > 0)
        {
            for (i = qn - 1; i >= 0; i--)
                q[i + limbs] = q[i];
            for (i = 0; i < limbs; i++)
                q[i] = 0;
            qn += limbs;
        }
        if (r > 0)
        {
            q[qn] = mpn_lshift(q, q, qn, r);
            qn += (q[qn] != 0);
        }
    }
    else if (sh < 0)
    {
        slong limbs = (-sh) / FLINT_BITS, r = (-sh) % FLINT_BITS;

        if (limbs >= qn)
        {
            *sticky = 1;
            return 0;
        }
        for (i = 0; i < limbs; i++)
            *sticky |= (q[i] != 0);
        if (limbs > 0)
        {
            for (i = 0; i < qn - limbs; i++)
                q[i] = q[i + limbs];
            qn -= limbs;
        }
        if (r > 0)
        {
            *sticky |= ((q[0] & ((UWORD(1) << r) - 1)) != 0);
            mpn_rshift(q, q, qn, r);
            qn -= (q[qn - 1] == 0);
        }
    }

    return qn;
}

/* Rounds (q + sticky epsilon) 10^-s (q nonzero, s a multiple of e) to
   res. Returns 1 on success and -1 on error. */
static int
_round_scaled(decfloat_t res, nn_srcptr q, slong qn, int negative, int sticky, slong s, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    nn_ptr rd = radix_integer_fit_limbs(&res->m, radix_set_mpn_need_alloc(qn, DECIMAL_CTX_RADIX(ctx)) + 1, DECIMAL_CTX_RADIX(ctx));
    slong rn = radix_set_mpn(rd, q, qn, DECIMAL_CTX_RADIX(ctx));
    fmpz sexp = -s / (slong) DECIMAL_CTX_E(ctx);
    return (_decfloat_set_round_limbs(res, rd, rn, negative, &sexp, sticky, prec, rnd, info, err, ctx) == GR_SUCCESS) ? 1 : -1;
}

static int
_decfloat_set_round_mpn_2exp_small(decfloat_t res, nn_srcptr mp, slong mn, int negative, slong t, slong s, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    ulong q[SMALL_LIMBS + 1];
    slong qn;
    int sticky;

    qn = _scaled_floor(q, &sticky, mp, mn, t, s);

    /* q = 0: bad exponent estimate */
    if (qn <= 0)
        return 0;

    return _round_scaled(res, q, qn, negative, sticky, s, prec, rnd, info, err, ctx);
}

/* Scaling for the conversion of x = m 2^t with p = bits(m) to prec +
   guard digits: sets s (a multiple of e) so that x 10^s has at least
   prec + guard digits (or is an integer), and returns the size in bits
   of the exact computation. */
static slong
_scaled_params(slong * sp, slong p, slong t, slong prec, slong guard, slong e)
{
    slong n2, E0, s;

    /* digits to compute: prec plus guard digits, so that the sticky bit
       falls well below the rounding position and the rounding error
       bound is tight unless the error is smaller than about
       10^-(prec+guard-1) relative to x */
    n2 = prec + guard;

    /* 10^E0 <= |x| < 10^(E0+2) approximately */
    E0 = (slong) floor((double) (p - 1 + t) * 0.30102999566398119521);

    /* scale by 10^s, s a multiple of e (s can be negative) */
    s = n2 - 1 - E0;
    if (s >= 0)
        s = e * ((s + e - 1) / e);
    else
        s = -e * ((-s) / e);

    /* x 10^s is already an integer for s >= -t, so a larger s (a huge
       requested precision) is pointless */
    if (s > 0)
    {
        slong s_exact = (t < 0) ? e * ((-t + e - 1) / e) : 0;
        if (s > s_exact)
            s = s_exact;
    }

    *sp = s;

    /* size of the exact computation */
    return p + (slong) (2.33 * FLINT_ABS(s)) + FLINT_MAX(t + s, 0);
}

/* Rounds x = m 2^t (m nonzero) to prec digits (prec != DECIMAL_PREC_EXACT),
   using exact integer arithmetic on a scaled truncation of x when the sizes
   are moderate. Returns 0 if the sizes are too large (the caller should
   then use arb). */
static int
_decfloat_set_round_fmpz_2exp_scaled_guard(decfloat_t res, const fmpz_t m, slong t, slong prec, int rnd,
    slong guard, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    slong p, s, sh, limit, size;
    fmpz_t q, r, u, am;
    fmpz sexp;
    int negative, sticky, status;

    limit = DECIMAL_EXACT_CONV_LIMIT(prec);
    p = fmpz_bits(m);
    negative = fmpz_sgn(m) < 0;

    if (FLINT_ABS(t) > limit)
        return 0;

    size = _scaled_params(&s, p, t, prec, guard, e);
    if (size > limit)
        return 0;

    if (size <= SMALL_LIMBS * FLINT_BITS - 2 && fmpz_size(m) <= 2)
    {
        int r;
        if (!COEFF_IS_MPZ(*m))
        {
            ulong m0 = FLINT_ABS(*m);
            r = _decfloat_set_round_mpn_2exp_small(res, &m0, 1, negative, t, s, prec, rnd, info, err, ctx);
        }
        else
        {
            __mpz_struct * mz = COEFF_TO_PTR(*m);
            r = _decfloat_set_round_mpn_2exp_small(res, mz->_mp_d, FLINT_ABS(mz->_mp_size), negative, t, s, prec, rnd, info, err, ctx);
        }
        if (r != 0)
            return r;
    }

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(u);
    fmpz_init(am);
    fmpz_abs(am, m);

    sticky = 0;

    if (s >= 0)
    {
        /* |x| 10^s = |m| 5^s 2^(t+s) */
        fmpz_ui_pow_ui(u, 5, s);
        fmpz_mul(q, u, am);
        sh = t + s;
        if (sh >= 0)
        {
            fmpz_mul_2exp(q, q, sh);
        }
        else
        {
            fmpz_fdiv_r_2exp(r, q, -sh);
            sticky = !fmpz_is_zero(r);
            fmpz_fdiv_q_2exp(q, q, -sh);
        }
    }
    else
    {
        /* |x| 10^s = |m| 2^(t+s) / 5^(-s) */
        fmpz_ui_pow_ui(u, 5, -s);
        sh = t + s;
        if (sh >= 0)
        {
            fmpz_mul_2exp(q, am, sh);
            fmpz_fdiv_qr(q, r, q, u);
            sticky = !fmpz_is_zero(r);
        }
        else
        {
            fmpz_fdiv_qr(q, r, am, u);
            sticky = !fmpz_is_zero(r);
            fmpz_fdiv_r_2exp(r, q, -sh);
            sticky |= !fmpz_is_zero(r);
            fmpz_fdiv_q_2exp(q, q, -sh);
        }
    }

    if (fmpz_is_zero(q))
    {
        /* cannot happen with the digit estimate above, but be safe */
        status = 0;
    }
    else
    {
        radix_integer_set_fmpz(&res->m, q, DECIMAL_CTX_RADIX(ctx));
        radix_integer_fit_limbs(&res->m, FLINT_ABS(res->m.size) + 1, DECIMAL_CTX_RADIX(ctx));
        sexp = -s / e;
        status = _decfloat_set_round_limbs(res, res->m.d, FLINT_ABS(res->m.size), negative, &sexp, sticky, prec, rnd, info, err, ctx);
        status = (status == GR_SUCCESS) ? 1 : -1;
    }

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(u);
    fmpz_clear(am);
    return status;
}

/* whether the error bound is not limited by the resolution of the
   truncation, i.e. err > 10^(E - prec - guard + 3) */
static int
_scaled_err_tight(const decfloat_t res, const decmag_t err, slong prec, slong guard, gr_ctx_t ctx)
{
    decmag_t u;
    int tight;
    _decmag_init(u, ctx);
    _decmag_set_10exp_si(u, _decfloat_sci_exp_clamped(res, ctx) - prec - guard + 3, ctx);
    tight = _decmag_cmp(err, u, ctx) > 0;
    _decmag_clear(u, ctx);
    return tight;
}

/* Rounds x = m 2^t (m nonzero) to prec digits (prec != DECIMAL_PREC_EXACT),
   using exact integer arithmetic on a scaled truncation of x when the sizes
   are moderate. Returns 0 if the sizes are too large (the caller should
   then use arb), -1 on error, 1 on success. When an error bound is
   requested and it is limited by the resolution of the truncation (the
   value is extremely close to a prec-digit number), the computation is
   repeated with more guard digits while the sizes allow. */
int
_decfloat_set_round_fmpz_2exp_scaled(decfloat_t res, const fmpz_t m, slong t, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    slong guard = 20, guard_max = 40 * prec + 4000;
    decimal_rounding_info info2;
    int r;

    if (info == NULL)
        info = &info2;

    r = _decfloat_set_round_fmpz_2exp_scaled_guard(res, m, t, prec, rnd, guard, info, err, ctx);

    if (r != 1 || err == NULL)
        return r;

    while (info->inexact && !DECMAG_IS_ZERO(err) && !DECFLOAT_IS_SPECIAL(res)
        && !info->underflow && !info->overflow && guard < guard_max
        && !_scaled_err_tight(res, err, prec, guard, ctx))
    {
        decfloat_t res3;
        decmag_t err3;
        decimal_rounding_info info3;

        guard *= 2;

        decfloat_init(res3, ctx);
        _decmag_init(err3, ctx);
        r = _decfloat_set_round_fmpz_2exp_scaled_guard(res3, m, t, prec, rnd, guard, &info3, err3, ctx);
        if (r == 1)
        {
            decfloat_swap(res, res3, ctx);
            _decmag_swap(err, err3, ctx);
            *info = info3;
        }
        decfloat_clear(res3, ctx);
        _decmag_clear(err3, ctx);

        /* size limit reached: keep the result we have */
        if (r != 1)
            break;
    }

    return 1;
}

int
_decfloat_set_round_fmpz_2exp_err(decfloat_t res, const fmpz_t m, const fmpz_t t, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    int status;

    if (fmpz_is_zero(m))
    {
        if (info != NULL)
            info->inexact = info->increased = info->underflow = info->overflow = 0;
        if (err != NULL)
            _decmag_zero(err, ctx);
        return decfloat_zero(res, ctx);
    }

    if (prec != DECIMAL_PREC_EXACT)
    {
        int r = fmpz_fits_si(t) ? _decfloat_set_round_fmpz_2exp_scaled(res, m, fmpz_get_si(t), prec, rnd, info, err, ctx) : 0;

        if (r == 1)
            return GR_SUCCESS;
        if (r == -1)
            return GR_UNABLE;
    }

    if (prec != DECIMAL_PREC_EXACT || fmpz_bits(t) > 20)
    {
        /* too large for an exact conversion: round via a scaled ball */
        arb_t a;

        if (prec == DECIMAL_PREC_EXACT)
            return GR_UNABLE;

        arb_init(a);
        arb_set_fmpz(a, m);
        arb_mul_2exp_fmpz(a, a, t);
        status = _decfloat_set_round_arb_err(res, err, a, prec, rnd, ctx);
        arb_clear(a);

        if (status == GR_SUCCESS && info != NULL)
        {
            info->inexact = (err != NULL) ? !DECMAG_IS_ZERO(err) : 1;
            info->increased = 0;   /* not tracked */
            info->underflow = 0;
            info->overflow = 0;
        }

        return status;
    }

    if (fmpz_sgn(t) >= 0)
    {
        fmpz_t u;
        fmpz_init(u);
        fmpz_mul_2exp(u, m, fmpz_get_ui(t));
        radix_integer_set_fmpz(&res->m, u, DECIMAL_CTX_RADIX(ctx));
        fmpz_clear(u);
        {
            fmpz zero = 0;
            slong n = FLINT_ABS(res->m.size);
            int negative = res->m.size < 0;
            radix_integer_fit_limbs(&res->m, n + 1, DECIMAL_CTX_RADIX(ctx));
            status = _decfloat_set_round_limbs(res, res->m.d, n, negative, &zero, 0, prec, rnd, info, err, ctx);
        }
    }
    else
    {
        /* m 2^-k = m 5^k 10^-k, with 10^-k = B^q 10^-rr */
        fmpz_t u, mm;
        slong e = DECIMAL_CTX_E(ctx);
        ulong k = -fmpz_get_si(t);
        slong q = -(slong) ((k + e - 1) / e);    /* floor(-k / e) */
        ulong rr = (ulong) (-q * e) - k;         /* -k - q e in [0, e) */
        fmpz qq;

        fmpz_init(u);
        fmpz_init(mm);
        fmpz_ui_pow_ui(u, 5, k);
        fmpz_mul(u, u, m);
        if (rr == 0)
            fmpz_swap(mm, u);
        else
            fmpz_mul_ui(mm, u, DECIMAL_CTX_RADIX(ctx)->bpow[rr]);
        radix_integer_set_fmpz(&res->m, mm, DECIMAL_CTX_RADIX(ctx));
        {
            slong n = FLINT_ABS(res->m.size);
            int negative = res->m.size < 0;
            radix_integer_fit_limbs(&res->m, n + 1, DECIMAL_CTX_RADIX(ctx));
            qq = q;
            status = _decfloat_set_round_limbs(res, res->m.d, n, negative, &qq, 0, prec, rnd, info, err, ctx);
        }
        fmpz_clear(u);
        fmpz_clear(mm);
    }

    return status;
}

int
decfloat_set_round_fmpz_2exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t t, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_set_round_fmpz_2exp_err(res, m, t, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_set_fmpz_2exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
{
    return decfloat_set_round_fmpz_2exp_fmpz(res, m, e, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_set_fmpz(decfloat_t res, const fmpz_t x, gr_ctx_t ctx)
{
    return decfloat_set_round_fmpz(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_set_ui(decfloat_t res, ulong x, gr_ctx_t ctx)
{
    return decfloat_set_round_ui(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_set_si(decfloat_t res, slong x, gr_ctx_t ctx)
{
    return decfloat_set_round_si(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

/* ------------------------------------------------------------------------- */
/*    From rationals, doubles, arf                                           */
/* ------------------------------------------------------------------------- */

int
decfloat_set_round_fmpq(decfloat_t res, const fmpq_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    decfloat_t a, b;
    int status;

    if (fmpz_is_one(fmpq_denref(x)))
        return decfloat_set_round_fmpz(res, fmpq_numref(x), prec, rnd, ctx);

    decfloat_init(a, ctx);
    decfloat_init(b, ctx);

    /* exact operands, which may be outside the exponent range */
    status = decfloat_set_round_fmpz(a, fmpq_numref(x), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
    status |= decfloat_set_round_fmpz(b, fmpq_denref(x), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);

    if (status == GR_SUCCESS)
        status = _decfloat_div(res, a, b, prec, rnd, NULL, NULL, ctx);

    decfloat_clear(a, ctx);
    decfloat_clear(b, ctx);

    return status;
}

int
decfloat_set_fmpq(decfloat_t res, const fmpq_t x, gr_ctx_t ctx)
{
    return decfloat_set_round_fmpq(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

/* Fast path for rounding a finite nonzero arf x with a short mantissa
   and a moderate exponent to prec != DECIMAL_PREC_EXACT digits (without
   an error bound). Returns 1 on success, -1 on error and 0 if not
   applicable. */
int
_decfloat_set_round_arf_small(decfloat_t res, const arf_t x, slong prec, int rnd, decimal_rounding_info * info, gr_ctx_t ctx)
{
    nn_srcptr xp;
    slong xn, t, s, size;

    if (ARF_SIZE(x) > SMALL_MN || arf_is_special(x) || COEFF_IS_MPZ(ARF_EXP(x))
        || FLINT_ABS(ARF_EXP(x)) >= (WORD(1) << 20))
        return 0;

    ARF_GET_MPN_READONLY(xp, xn, x);
    t = ARF_EXP(x) - xn * FLINT_BITS;
    /* no error bound is needed, so two guard digits suffice */
    size = _scaled_params(&s, xn * FLINT_BITS, t, prec, 2, DECIMAL_CTX_E(ctx));

    if (size > SMALL_LIMBS * FLINT_BITS - 2)
        return 0;

    return _decfloat_set_round_mpn_2exp_small(res, xp, xn, ARF_SGNBIT(x), t, s, prec, rnd, info, NULL, ctx);
}

/* Rounds the ball [mid +/- rad] (mid nonzero, rad nonzero, both with
   small exponents and a short mid) to prec digits using a
   single scaled integer conversion: with Q = floor(|mid| 10^s) and
   R = ceil(rad 10^s), |x| 10^s lies in [Q - R, Q + R + 1], where s is
   chosen as described below. Returns 1
   if both ends round to the same value (set in res), 0 if undecided
   and -2 if not applicable. */
int
_decfloat_round_arb_small(decfloat_t res, const arf_t mid, const mag_t rad, slong prec, int rnd, gr_ctx_t ctx)
{
    ulong Q[SMALL_LIMBS + 2], R[SMALL_LIMBS + 2], m;
    nn_srcptr xp;
    slong xn, t, s, size, qn, rn, d, guard;
    int sq, sr, r1, r2, negative;
    decfloat_t rh;

    if (ARF_SIZE(mid) > SMALL_MN || COEFF_IS_MPZ(ARF_EXP(mid)) || COEFF_IS_MPZ(MAG_EXP(rad))
        || FLINT_ABS(ARF_EXP(mid)) >= (WORD(1) << 20) || FLINT_ABS(MAG_EXP(rad)) >= (WORD(1) << 20))
        return -2;

    d = ARF_EXP(mid) - MAG_EXP(rad);
    if (d < 8)
        return -2;

    ARF_GET_MPN_READONLY(xp, xn, mid);
    negative = ARF_SGNBIT(mid);
    t = ARF_EXP(mid) - xn * FLINT_BITS;

    /* Scale by 10^s (s a multiple of e) so that the unit is comparable
       to rad (about d log10(2) digits relative to mid; this resolution
       is not limited by the exactness of mid), rounding s down to a
       multiple of e when that leaves at least prec + 8 digits: the
       widening by a few units only slightly increases the chance of an
       undecided result, while the sizes stay small. */
    {
        slong e = DECIMAL_CTX_E(ctx), E0, st;

        guard = FLINT_MAX((slong) (d * 0.30102999566398119521) - prec, 2);
        /* 10^E0 <= |mid| < 10^(E0+2) */
        E0 = (slong) floor((double) (ARF_EXP(mid) - 1) * 0.30102999566398119521);
        st = prec + guard - 1 - E0;
        s = (st >= 0) ? e * (st / e) : -e * ((-st + e - 1) / e);
        if (E0 + s + 1 < prec + 8)
            s += e;
        size = xn * FLINT_BITS + (slong) (2.33 * FLINT_ABS(s)) + FLINT_MAX(t + s, 0);
        if (s < 0 || size > SMALL_LIMBS * FLINT_BITS - 2)
            return -2;
    }

    qn = _scaled_floor(Q, &sq, xp, xn, t, s);
    if (qn <= 0)
        return -2;

    m = MAG_MAN(rad);
    rn = _scaled_floor(R, &sr, &m, 1, MAG_EXP(rad) - MAG_BITS, s);
    if (rn < 0 || rn > qn)
        return (rn < 0) ? -2 : 0;

    /* R = ceil(rad 10^s) + 1 when mid 10^s was inexact, so that
       C = Q + R bounds |x| 10^s from above */
    R[rn] = 0;
    if (sr + sq != 0)
    {
        if (rn == 0)
        {
            R[0] = sr + sq;
            rn = 1;
        }
        else
        {
            R[rn] = mpn_add_1(R, R, rn, sr + sq);
            rn += (R[rn] != 0);
        }
    }

    /* A = Q - ceil(rad 10^s) <= |x| 10^s */
    {
        ulong A[SMALL_LIMBS + 2], C[SMALL_LIMBS + 2];
        slong an, cn;

        flint_mpn_copyi(C, Q, qn);
        C[qn] = (rn == 0) ? 0 : mpn_add(C, C, qn, R, rn);
        cn = qn + (C[qn] != 0);

        /* undo the inexactness unit of mid for the lower end */
        if (sq != 0)
        {
            if (rn == 0 || mpn_sub_1(R, R, rn, 1))
                return -2;
            while (rn > 0 && R[rn - 1] == 0)
                rn--;
        }

        if (rn > qn || (rn == qn && mpn_cmp(Q, R, qn) <= 0))
            return 0;   /* the interval reaches zero */

        flint_mpn_copyi(A, Q, qn);
        if (rn != 0)
            mpn_sub(A, A, qn, R, rn);
        an = qn;
        while (A[an - 1] == 0)
            an--;

        r1 = _round_scaled(res, A, an, negative, 0, s, prec, rnd, NULL, NULL, ctx);
        if (r1 != 1)
            return 0;

        decfloat_init(rh, ctx);
        r2 = _round_scaled(rh, C, cn, negative, 0, s, prec, rnd, NULL, NULL, ctx);
        r1 = (r2 == 1) && (decfloat_equal(res, rh, ctx) == T_TRUE);
        decfloat_clear(rh, ctx);
    }

    return r1;
}

int
decfloat_set_round_arf(decfloat_t res, const arf_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz_t m, e;
    int status;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return decfloat_zero(res, ctx);
        if (arf_is_pos_inf(x))
            return decfloat_pos_inf(res, ctx);
        if (arf_is_neg_inf(x))
            return decfloat_neg_inf(res, ctx);
        return decfloat_nan(res, ctx);
    }

    if (prec != DECIMAL_PREC_EXACT)
    {
        int r = _decfloat_set_round_arf_small(res, x, prec, rnd, NULL, ctx);
        if (r != 0)
            return (r == 1) ? GR_SUCCESS : GR_UNABLE;
    }

    fmpz_init(m);
    fmpz_init(e);
    arf_get_fmpz_2exp(m, e, x);
    status = decfloat_set_round_fmpz_2exp_fmpz(res, m, e, prec, rnd, ctx);
    fmpz_clear(m);
    fmpz_clear(e);
    return status;
}

int
decfloat_set_arf(decfloat_t res, const arf_t x, gr_ctx_t ctx)
{
    return decfloat_set_round_arf(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_set_round_d(decfloat_t res, double x, slong prec, int rnd, gr_ctx_t ctx)
{
    arf_t t;
    int status;
    arf_init(t);
    arf_set_d(t, x);
    status = decfloat_set_round_arf(res, t, prec, rnd, ctx);
    arf_clear(t);
    return status;
}

int
decfloat_set_d(decfloat_t res, double x, gr_ctx_t ctx)
{
    return decfloat_set_round_d(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

/* ------------------------------------------------------------------------- */
/*    To integers and rationals                                              */
/* ------------------------------------------------------------------------- */

int
decfloat_get_fmpz_10exp_fmpz(fmpz_t m, fmpz_t t, const decfloat_t x, gr_ctx_t ctx)
{
    ulong v0;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
        {
            fmpz_zero(m);
            fmpz_zero(t);
            return GR_SUCCESS;
        }
        return GR_DOMAIN;
    }

    radix_integer_get_fmpz(m, &x->m, DECIMAL_CTX_RADIX(ctx));
    fmpz_mul_ui(t, &x->exp, DECIMAL_CTX_E(ctx));

    v0 = _radix_valuation_digits_1(x->m.d[0], DECIMAL_CTX_RADIX(ctx));
    if (v0 != 0)
    {
        fmpz_divexact_ui(m, m, DECIMAL_CTX_RADIX(ctx)->bpow[v0]);
        fmpz_add_ui(t, t, v0);
    }

    return GR_SUCCESS;
}

int
_decfloat_set_scalar_exact(decfloat_t res, const void * y, int type, gr_ctx_t ctx)
{
    int rnd = DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS;

    switch (type)
    {
        case DECIMAL_SCALAR_UI: return decfloat_set_round_ui(res, *(const ulong *) y, DECIMAL_PREC_EXACT, rnd, ctx);
        case DECIMAL_SCALAR_SI: return decfloat_set_round_si(res, *(const slong *) y, DECIMAL_PREC_EXACT, rnd, ctx);
        case DECIMAL_SCALAR_FMPZ: return decfloat_set_round_fmpz(res, (const fmpz *) y, DECIMAL_PREC_EXACT, rnd, ctx);
        default: return decfloat_set_round_d(res, *(const double *) y, DECIMAL_PREC_EXACT, rnd, ctx);
    }
}

int
decfloat_get_fmpz(fmpz_t res, const decfloat_t x, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
        {
            fmpz_zero(res);
            return GR_SUCCESS;
        }
        return GR_DOMAIN;
    }

    if (fmpz_sgn(&x->exp) < 0)
        return GR_DOMAIN;

    if (COEFF_IS_MPZ(x->exp) || x->exp > (slong) (FLINT_MIN(ctx->size_limit, DECIMAL_CONV_DIGITS_LIMIT) / e))
        return GR_UNABLE;

    radix_integer_get_fmpz(res, &x->m, DECIMAL_CTX_RADIX(ctx));

    if (!fmpz_is_zero(&x->exp))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_ui_pow_ui(t, 10, e * x->exp);
        fmpz_mul(res, res, t);
        fmpz_clear(t);
    }

    return GR_SUCCESS;
}

int
decfloat_get_fmpq(fmpq_t res, const decfloat_t x, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
        {
            fmpq_zero(res);
            return GR_SUCCESS;
        }
        return GR_DOMAIN;
    }

    if (fmpz_sgn(&x->exp) >= 0)
    {
        int status = decfloat_get_fmpz(fmpq_numref(res), x, ctx);
        fmpz_one(fmpq_denref(res));
        return status;
    }

    if (COEFF_IS_MPZ(x->exp) || -x->exp > (slong) (FLINT_MIN(ctx->size_limit, DECIMAL_CONV_DIGITS_LIMIT) / e))
        return GR_UNABLE;

    radix_integer_get_fmpz(fmpq_numref(res), &x->m, DECIMAL_CTX_RADIX(ctx));
    fmpz_ui_pow_ui(fmpq_denref(res), 10, e * (-x->exp));
    fmpq_canonicalise(res);
    return GR_SUCCESS;
}

int
decfloat_get_si(slong * res, const decfloat_t x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init(t);
    status = decfloat_get_fmpz(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_fits_si(t))
            *res = fmpz_get_si(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

int
decfloat_get_ui(ulong * res, const decfloat_t x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init(t);
    status = decfloat_get_fmpz(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(t) >= 0 && fmpz_abs_fits_ui(t))
            *res = fmpz_get_ui(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

/* ------------------------------------------------------------------------- */
/*    To arf and double                                                      */
/* ------------------------------------------------------------------------- */

/* exact conversion of m * 10^t (t < 0) if 5^(-t) divides m; returns 1 on
   success */
static int
_get_arf_exact_neg(arf_t res, const fmpz_t m, slong t, slong prec_bits, int rnd)
{
    fmpz_t q, r, e5;
    int ok;

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(e5);
    fmpz_ui_pow_ui(e5, 5, -t);
    fmpz_fdiv_qr(q, r, m, e5);
    ok = fmpz_is_zero(r);
    if (ok)
    {
        fmpz_set_si(e5, t);
        arf_set_fmpz_2exp(res, q, e5);
        if (prec_bits != ARF_PREC_EXACT)
            arf_set_round(res, res, prec_bits, rnd);
    }
    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(e5);
    return ok;
}

/* A cheap necessary condition for a finite x = M 10^(e v) to be a dyadic
   rational: for v < 0, 5^(-e v) must divide M, and M = d_0 mod 5^e. */
int
_decfloat_maybe_dyadic(const decfloat_t x, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    ulong d, p;
    slong j, k;

    if (DECFLOAT_IS_SPECIAL(x) || fmpz_sgn(&x->exp) >= 0)
        return 1;

    k = COEFF_IS_MPZ(x->exp) ? e : FLINT_MIN(e, -x->exp * e);
    d = x->m.d[0];
    if (k <= DECIMAL_POW5_MAX_EXP)
        return (d % _pow5_tab[k]) == 0;
    for (j = 0, p = 1; j < k; j++)
        p *= 5;
    return (d % p) == 0;
}

int
decfloat_get_arf(arf_t res, const decfloat_t x, slong prec_bits, int rnd, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    fmpz_t m, t;
    int status = GR_SUCCESS;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
            arf_zero(res);
        else if (DECFLOAT_IS_POS_INF(x))
            arf_pos_inf(res);
        else if (DECFLOAT_IS_NEG_INF(x))
            arf_neg_inf(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }

    /* not a dyadic rational: avoid the exact conversion attempt */
    if (prec_bits == ARF_PREC_EXACT && !_decfloat_maybe_dyadic(x, ctx))
        return GR_UNABLE;

    fmpz_init(m);
    fmpz_init(t);

    radix_integer_get_fmpz(m, &x->m, DECIMAL_CTX_RADIX(ctx));
    fmpz_mul_ui(t, &x->exp, e);

    if (fmpz_is_zero(t))
    {
        if (prec_bits == ARF_PREC_EXACT)
            arf_set_fmpz(res, m);
        else
            arf_set_round_fmpz(res, m, prec_bits, rnd);
    }
    else if (fmpz_bits(t) <= 16)
    {
        slong tt = fmpz_get_si(t);
        fmpz_t p;
        fmpz_init(p);

        if (tt > 0)
        {
            fmpz_ui_pow_ui(p, 10, tt);
            fmpz_mul(m, m, p);
            if (prec_bits == ARF_PREC_EXACT)
                arf_set_fmpz(res, m);
            else
                arf_set_round_fmpz(res, m, prec_bits, rnd);
        }
        else if (prec_bits == ARF_PREC_EXACT)
        {
            if (!_get_arf_exact_neg(res, m, tt, prec_bits, rnd))
                status = GR_UNABLE;
        }
        else
        {
            arf_t a, b;
            fmpz_ui_pow_ui(p, 10, -tt);
            arf_init(a);
            arf_init(b);
            arf_set_fmpz(a, m);
            arf_set_fmpz(b, p);
            arf_div(res, a, b, prec_bits, rnd);
            arf_clear(a);
            arf_clear(b);
        }

        fmpz_clear(p);
    }
    else
    {
        /*
            Huge decimal exponent. A value that is exactly representable at
            the target precision (or lies on a rounding boundary) cannot be
            handled by Ziv's strategy, so detect that case first: for t > 0
            the value m 10^t has about bits(m) + 2.32 t bits, for t < 0 it is
            representable iff 5^(-t) divides m, which requires m >= 5^(-t).
        */
        double tbits = fmpz_get_d(t) * 2.3219280948873623479;
        slong mbits = fmpz_bits(m);
        int decided = 0;

        if (fmpz_sgn(t) > 0)
        {
            if (prec_bits == ARF_PREC_EXACT || tbits + mbits <= prec_bits + 2)
            {
                if (tbits + mbits > (double) DECIMAL_CONV_DIGITS_LIMIT * 3.33)
                {
                    status = GR_UNABLE;
                }
                else
                {
                    fmpz_t p;
                    fmpz_init(p);
                    fmpz_ui_pow_ui(p, 10, fmpz_get_ui(t));
                    fmpz_mul(m, m, p);
                    if (prec_bits == ARF_PREC_EXACT)
                        arf_set_fmpz(res, m);
                    else
                        arf_set_round_fmpz(res, m, prec_bits, rnd);
                    fmpz_clear(p);
                }
                decided = 1;
            }
        }
        else
        {
            if (mbits + 2 >= -tbits)
            {
                /* divisibility by 5^(-t) is possible: check exactly */
                if (_get_arf_exact_neg(res, m, fmpz_get_si(t), prec_bits, rnd))
                    decided = 1;
                else if (prec_bits == ARF_PREC_EXACT)
                {
                    status = GR_UNABLE;
                    decided = 1;
                }
            }
            else if (prec_bits == ARF_PREC_EXACT)
            {
                status = GR_UNABLE;
                decided = 1;
            }
        }

        if (!decided)
        {
            /* not exactly representable (nor a tie): Ziv's strategy */
            arb_t a, b;
            slong wp;

            arb_init(a);
            arb_init(b);

            for (wp = prec_bits + 32; ; wp *= 2)
            {
                arb_set_fmpz(a, m);
                _decimal_arb_10_pow_fmpz(b, t, wp);
                arb_mul(a, a, b, wp);

                if (arb_can_round_arf(a, prec_bits, rnd))
                {
                    arf_set_round(res, arb_midref(a), prec_bits, rnd);
                    break;
                }

                /* the distance to a rounding boundary is at least about
                   2^-(mbits + prec_bits) relative to the value */
                if (wp > 100 * prec_bits + 1000000 && wp > 2 * (mbits + prec_bits) + 128)
                {
                    status = GR_UNABLE;
                    break;
                }
            }

            arb_clear(a);
            arb_clear(b);
        }
    }

    fmpz_clear(m);
    fmpz_clear(t);
    return status;
}

int
decfloat_get_d(double * res, const decfloat_t x, gr_ctx_t ctx)
{
    arf_t t;
    slong E;
    int status;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
            *res = 0.0;
        else if (DECFLOAT_IS_POS_INF(x))
            *res = D_INF;
        else if (DECFLOAT_IS_NEG_INF(x))
            *res = -D_INF;
        else
            *res = D_NAN;
        return GR_SUCCESS;
    }

    /* out of the normal double range: refuse rather than returning
       0, a subnormal, or inf */
    if (!decfloat_get_sci_exp_si(&E, x, ctx) || E > 308 || E < -307)
        return GR_UNABLE;

    arf_init(t);
    status = decfloat_get_arf(t, x, 53, ARF_RND_NEAR, ctx);
    *res = arf_get_d(t, ARF_RND_NEAR);
    arf_clear(t);
    return status;
}

int
decfloat_get_fmpz_fixed_si(fmpz_t res, const decfloat_t x, slong e, int rnd, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;

    decfloat_init(t, ctx);
    status = decfloat_mul_10exp_si_round(t, x, -e, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_round_to_int(t, t, rnd, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, NULL, NULL, ctx);
    if (status == GR_SUCCESS)
        status = decfloat_get_fmpz(res, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

/* ------------------------------------------------------------------------- */
/*    Generic conversion                                                     */
/* ------------------------------------------------------------------------- */

/* Conversion between decimal contexts with possibly different limb sizes. */
int
_decfloat_set_decfloat_other(decfloat_t res, const decfloat_t y, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    fmpz_t m, t;
    int status;

    if (DECIMAL_CTX_E(x_ctx) == DECIMAL_CTX_E(ctx))
        return decfloat_set(res, y, ctx);

    if (DECFLOAT_IS_SPECIAL(y))
    {
        if (DECFLOAT_IS_ZERO(y)) return decfloat_zero(res, ctx);
        if (DECFLOAT_IS_POS_INF(y)) return decfloat_pos_inf(res, ctx);
        if (DECFLOAT_IS_NEG_INF(y)) return decfloat_neg_inf(res, ctx);
        return decfloat_nan(res, ctx);
    }

    fmpz_init(m);
    fmpz_init(t);
    status = decfloat_get_fmpz_10exp_fmpz(m, t, y, x_ctx);
    if (status == GR_SUCCESS)
        status = decfloat_set_fmpz_10exp_fmpz(res, m, t, ctx);
    fmpz_clear(m);
    fmpz_clear(t);
    return status;
}

int
decfloat_set_other(decfloat_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return decfloat_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return decfloat_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
            return decfloat_set_arf(res, x, ctx);

        case GR_CTX_RR_ARB:
            if (DECIMAL_CTX_IS_EXACT(ctx) && !arb_is_exact((arb_srcptr) x))
                return GR_UNABLE;
            return decfloat_set_arf(res, arb_midref((arb_srcptr) x), ctx);

        case GR_CTX_DECFLOAT:
            return _decfloat_set_decfloat_other(res, x, x_ctx, ctx);

        case GR_CTX_DECBALL:
            if (DECIMAL_CTX_IS_EXACT(ctx) && !DECMAG_IS_ZERO(DECBALL_RADREF((decball_srcptr) x)))
                return GR_UNABLE;
            return _decfloat_set_decfloat_other(res, DECBALL_MIDREF((decball_srcptr) x), x_ctx, ctx);

        case GR_CTX_DECCFLOAT:
            if (!_deccfloat_is_real((deccfloat_srcptr) x))
                return _deccfloat_is_nan((deccfloat_srcptr) x) ? GR_UNABLE : GR_DOMAIN;
            return _decfloat_set_decfloat_other(res, DECCFLOAT_REALREF((deccfloat_srcptr) x), x_ctx, ctx);

        case GR_CTX_DECCBALL:
            if (!_deccball_is_real((deccball_srcptr) x, x_ctx))
                return _decball_contains_zero(DECCBALL_IMAGREF((deccball_srcptr) x), x_ctx) ? GR_UNABLE : GR_DOMAIN;
            if (DECIMAL_CTX_IS_EXACT(ctx) && !DECMAG_IS_ZERO(DECBALL_RADREF(DECCBALL_REALREF((deccball_srcptr) x))))
                return GR_UNABLE;
            return _decfloat_set_decfloat_other(res, DECBALL_MIDREF(DECCBALL_REALREF((deccball_srcptr) x)), x_ctx, ctx);

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            return decfloat_set_qqbar(res, x, ctx);

        default:
            {
                gr_ctx_t cctx;
                acb_t z;
                int status;

                gr_ctx_init_complex_acb(cctx, 20 + (DECIMAL_CTX_IS_EXACT(ctx) ? 64 : 4 * DECIMAL_CTX_PREC(ctx)));
                acb_init(z);

                status = gr_set_other(z, x, x_ctx, cctx);

                if (status == GR_SUCCESS)
                {
                    if (!acb_is_real(z))
                        status = GR_DOMAIN;
                    else if (DECIMAL_CTX_IS_EXACT(ctx) && !acb_is_exact(z))
                        status = GR_UNABLE;
                    else
                        status = decfloat_set_arf(res, arb_midref(acb_realref(z)), ctx);
                }

                acb_clear(z);
                gr_ctx_clear(cctx);

                return status;
            }
    }
}
