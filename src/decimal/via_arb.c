/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "decimal.h"
#include "mag.h"
#include "fmpz_extras.h"
#include "fmpq.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"

/* Not performance-critical: optimize for size. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* bits needed for prec decimal digits, with guard bits */
slong
_decimal_digits_to_bits(slong prec)
{
    return (slong) (prec * 3.3219280948873623479) + 30;
}

/* ------------------------------------------------------------------------- */
/*    Conversions                                                            */
/* ------------------------------------------------------------------------- */

/*
    Rigorous enclosure of a decimal number in an arb ball at prec_bits bits.

    Small operands (whose exact binary conversion costs no more than the
    target precision) are converted exactly when they are dyadic rationals,
    so that exact inputs give exact balls. Otherwise only the leading limbs
    of the mantissa are used: with k limbs kept, x = (M + f) B^s where M is
    the integer formed by the top k limbs and 0 <= f < 1 (f > 0 iff limbs
    were dropped, since the lowest limb of a canonical mantissa is nonzero),
    so x lies in [M +/- 1] B^s. The cost is then independent of the length
    of the mantissa.
*/
int
decfloat_get_arb(arb_t res, const decfloat_t x, slong prec_bits, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    slong n, k, T;
    double size_bits;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
            arb_zero(res);
        else if (DECFLOAT_IS_POS_INF(x))
            arb_pos_inf(res);
        else if (DECFLOAT_IS_NEG_INF(x))
            arb_neg_inf(res);
        else
            arb_indeterminate(res);
        return GR_SUCCESS;
    }

    n = FLINT_ABS(x->m.size);

    /* size of the exact binary representation of M 10^t, in bits */
    if (COEFF_IS_MPZ(x->exp) || FLINT_ABS(x->exp) > WORD_MAX / (4 * e) - n)
        T = WORD_MAX / 4;
    else
        T = FLINT_ABS(x->exp) * e;
    size_bits = ((double) n * e + (double) T) * 3.3219280948873623479;

    if (size_bits <= 4.0 * (double) prec_bits + 4096.0)
    {
        if (_decfloat_maybe_dyadic(x, ctx) && decfloat_get_arf(arb_midref(res), x, ARF_PREC_EXACT, ARF_RND_DOWN, ctx) == GR_SUCCESS)
        {
            /* exact when the value is a dyadic rational */
            mag_zero(arb_radref(res));
            arb_set_round(res, res, prec_bits);
        }
        else
        {
            /* correctly rounded toward zero: the error is less than one ulp */
            fmpz_t t;
            GR_MUST_SUCCEED(decfloat_get_arf(arb_midref(res), x, prec_bits, ARF_RND_DOWN, ctx));
            fmpz_init(t);
            fmpz_sub_ui(t, ARF_EXPREF(arb_midref(res)), prec_bits);
            mag_one(arb_radref(res));
            mag_mul_2exp_fmpz(arb_radref(res), arb_radref(res), t);
            fmpz_clear(t);
        }
    }
    else
    {
        radix_integer_struct top;
        fmpz_t M, s;
        arb_t p;

        /* k limbs carry at least prec_bits + 64 bits */
        k = (slong) ((prec_bits + 64) / (e * 3.3219280948873623479)) + 2;
        k = FLINT_MIN(k, n);

        fmpz_init(M);
        fmpz_init(s);
        arb_init(p);

        top.d = x->m.d + (n - k);
        top.size = k;
        top.alloc = k;
        radix_integer_get_fmpz(M, &top, DECIMAL_CTX_RADIX(ctx));

        arf_set_fmpz(arb_midref(res), M);
        if (k < n)
            mag_one(arb_radref(res));
        else
            mag_zero(arb_radref(res));

        /* scale by 10^(e (exp + n - k)) */
        fmpz_add_ui(s, &x->exp, n - k);
        fmpz_mul_ui(s, s, e);
        _decimal_arb_10_pow_fmpz(p, s, prec_bits + 20);
        arb_mul(res, res, p, prec_bits);

        if (x->m.size < 0)
            arb_neg(res, res);

        fmpz_clear(M);
        fmpz_clear(s);
        arb_clear(p);
    }

    return GR_SUCCESS;
}

int
decball_get_arb(arb_t res, const decball_t x, slong prec_bits, gr_ctx_t ctx)
{
    decfloat_get_arb(res, &x->mid, prec_bits, ctx);

    if (!DECMAG_IS_ZERO(&x->rad))
    {
        mag_t r;
        mag_init(r);
        _decmag_get_mag(r, &x->rad, ctx);
        arb_add_error_mag(res, r);
        mag_clear(r);
    }

    return GR_SUCCESS;
}

/* Rigorous conversion of an arb ball with a moderate exponent to a decimal
   ball at the context precision, by exact conversion of the midpoint. */
static int
_decball_set_arb_direct(decball_t res, const arb_t x, gr_ctx_t ctx)
{
    decimal_rounding_info info;
    decmag_t err;
    int status;

    _decmag_init(err, ctx);

    if (DECIMAL_CTX_PREC(ctx) == DECIMAL_PREC_EXACT)
    {
        status = decfloat_set_round_arf(&res->mid, arb_midref(x), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx);
        info.inexact = info.increased = info.underflow = info.overflow = 0;
        _decmag_zero(err, ctx);
    }
    else if (arf_is_special(arb_midref(x)))
    {
        status = decfloat_set_round_arf(&res->mid, arb_midref(x), DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
        info.inexact = info.increased = info.underflow = info.overflow = 0;
        _decmag_zero(err, ctx);
    }
    else
    {
        fmpz_t m, t;
        fmpz_init(m);
        fmpz_init(t);
        arf_get_fmpz_2exp(m, t, arb_midref(x));
        status = _decfloat_set_round_fmpz_2exp_err(&res->mid, m, t, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx),
            &info, DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx);
        fmpz_clear(m);
        fmpz_clear(t);
    }

    if (status == GR_SUCCESS)
    {
        _decmag_set_mag(&res->rad, arb_radref(x), ctx);
        _decball_add_rounding_error(res, &info, err, DECIMAL_CTX_PREC(ctx), ctx);
    }

    _decmag_clear(err, ctx);
    return status;
}

/* E approximately floor(t log10(2)) for a binary exponent t, computed
   accurately enough that the residual exponent is small */
static void
_approx_decimal_exponent(fmpz_t E, const fmpz_t t)
{
    arb_t u, v;
    slong prec = fmpz_bits(t) + 30;

    if (!COEFF_IS_MPZ(*t))
    {
        /* fixed-point product with floor(log10(2) 2^FLINT_BITS) */
        slong tt = *t;
        ulong hi, lo;
#if FLINT_BITS == 64
        umul_ppmm(hi, lo, (ulong) FLINT_ABS(tt), UWORD(0x4d104d427de7fbcc));
#else
        umul_ppmm(hi, lo, (ulong) FLINT_ABS(tt), UWORD(0x4d104d42));
#endif
        (void) lo;
        if (tt >= 0)
            fmpz_set_ui(E, hi);
        else
            fmpz_neg_uiui(E, 0, hi + 1);
        return;
    }

    arb_init(u);
    arb_init(v);
    arb_const_log2(u, prec);
    arb_const_log10(v, prec);
    arb_div(u, u, v, prec);
    arb_mul_fmpz(u, u, t, prec);
    arf_get_fmpz(E, arb_midref(u), ARF_RND_FLOOR);
    arb_clear(u);
    arb_clear(v);
}

/* res = 10^t to wp bits, in time essentially linear in the size of t
   (arb_pow_fmpz is slow for huge t) */
void
_decimal_arb_10_pow_fmpz(arb_t res, const fmpz_t t, slong wp)
{
    if (fmpz_bits(t) <= 64)
    {
        arb_set_ui(res, 10);
        arb_pow_fmpz(res, res, t, wp);
    }
    else
    {
        /* 10^t = 2^N 2^f, u = t log2(10) = N + f */
        arb_t u, v;
        arf_t w;
        fmpz_t N;
        slong P = fmpz_bits(t) + wp + 20;

        arb_init(u);
        arb_init(v);
        arf_init(w);
        fmpz_init(N);

        arb_const_log10(u, P);
        arb_const_log2(v, P);
        arb_div(u, u, v, P);
        arb_mul_fmpz(u, u, t, P);
        arb_get_lbound_arf(w, u, P);
        arf_get_fmpz(N, w, ARF_RND_FLOOR);
        arb_sub_fmpz(u, u, N, wp + 20);
        /* 2^f = exp(f log 2) */
        arb_const_log2(v, wp + 20);
        arb_mul(u, u, v, wp + 20);
        arb_exp(res, u, wp + 10);
        arb_mul_2exp_fmpz(res, res, N);

        arb_clear(u);
        arb_clear(v);
        arf_clear(w);
        fmpz_clear(N);
    }
}

/* Conversion of a ball with a huge exponent, by scaling with a power of
   ten computed to wp bits. If eps is not NULL, it is set to a bound for
   the contribution of the scaling error to the radius. */
static int
_decball_set_arb_scaled(decball_t res, const arb_t x, slong wp, decmag_ptr eps, gr_ctx_t ctx)
{
    fmpz_t E;
    arb_t y, p;
    int status;

    fmpz_init(E);
    arb_init(y);
    arb_init(p);

    _approx_decimal_exponent(E, ARF_EXPREF(arb_midref(x)));

    /* y = x 10^-E, as a product or a quotient by 10^|E| */
    if (fmpz_abs_fits_ui(E))
    {
        fmpz_t absE;
        fmpz_init(absE);
        fmpz_abs(absE, E);
        arb_ui_pow_ui(p, 10, fmpz_get_ui(absE), wp);
        if (fmpz_sgn(E) <= 0)
            arb_mul(y, x, p, wp);
        else
            arb_div(y, x, p, wp);
        fmpz_clear(absE);
    }
    else
    {
        fmpz_neg(E, E);
        _decimal_arb_10_pow_fmpz(p, E, wp);
        fmpz_neg(E, E);
        arb_mul(y, x, p, wp);
    }

    if (eps != NULL)
    {
        /* the relative error of the scaling is at most rad(p)/|p| plus
           the rounding error 2^(1-wp): eps = |y| (rad(p)/|p| + 2^(1-wp)) */
        mag_t t, u;
        mag_init(t);
        mag_init(u);
        arb_get_mag_lower(u, p);
        mag_div(t, arb_radref(p), u);
        mag_one(u);
        mag_mul_2exp_si(u, u, 1 - wp);
        mag_add(t, t, u);
        arb_get_mag(u, y);
        mag_mul(t, t, u);
        _decmag_set_mag(eps, t, ctx);
        mag_clear(t);
        mag_clear(u);
    }

    status = _decball_set_arb_direct(res, y, ctx);
    if (status == GR_SUCCESS)
        status = decball_mul_10exp_fmpz(res, res, E, ctx);
    if (status == GR_SUCCESS && eps != NULL)
        _decmag_mul_10exp_fmpz(eps, eps, E, ctx);

    fmpz_clear(E);
    arb_clear(y);
    arb_clear(p);
    return status;
}

/* If x = m 2^t (m odd, t >= 0, x large) is exactly representable with at
   most prec digits, sets res to x and returns 1. Otherwise returns 0. */
int
_decfloat_set_exact_decimal_fmpz_2exp(decfloat_t res, const fmpz_t m, const fmpz_t t, slong prec, gr_ctx_t ctx)
{
    slong p, tt, E0, D;
    fmpz_t u, q;
    int r = 0;

    if (fmpz_sgn(t) < 0 || !fmpz_fits_si(t) || prec == DECIMAL_PREC_EXACT)
        return 0;

    tt = fmpz_get_si(t);
    p = fmpz_bits(m);

    /* 10^E0 <= |x| approximately; x = D' 10^D requires 5^D | m and D <= t */
    E0 = (slong) floor((double) (p - 1 + tt) * 0.30102999566398119521) - 1;
    D = E0 + 1 - prec;

    if (D <= 0 || D > tt)
        return 0;

    /* 5^D has more than p bits */
    if ((double) D * 2.3219280948873623 > p + 1)
        return 0;

    /* cheap necessary condition before forming 5^D */
    if (D >= DECIMAL_POW5_MAX_EXP && !fmpz_divisible_ui(m, DECIMAL_POW5_MAX))
        return 0;

    fmpz_init(u);
    fmpz_init(q);
    fmpz_ui_pow_ui(u, 5, D);

    if (fmpz_divisible(m, u))
    {
        fmpz_t s;
        fmpz_init(s);
        fmpz_divexact(q, m, u);
        fmpz_set_si(s, tt - D);
        /* q 2^(t-D) has about prec + 2 digits */
        if (decfloat_set_round_fmpz_2exp_fmpz(res, q, s, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx) == GR_SUCCESS
            && decfloat_digits(res, ctx) <= prec)
        {
            fmpz_set_si(s, D);
            r = (decfloat_mul_10exp_fmpz_round(res, res, s, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx) == GR_SUCCESS);
        }
        fmpz_clear(s);
    }

    fmpz_clear(u);
    fmpz_clear(q);
    return r;
}

int
decball_set_arb(decball_t res, const arb_t x, gr_ctx_t ctx)
{
    fmpz_t m, t;
    decimal_rounding_info info;
    decmag_t err, eps;
    slong prec = DECIMAL_CTX_PREC(ctx);
    slong wp, wp_max;
    int r, status;

    /* balls represent real numbers: an infinite midpoint is not
       admitted (an infinite radius is, denoting the whole real line) */
    if (arf_is_nan(arb_midref(x)))
        return GR_UNABLE;
    if (arf_is_inf(arb_midref(x)))
        return GR_DOMAIN;
    if (mag_is_inf(arb_radref(x)))
        return decball_zero_pm_inf(res, ctx);

    if (arf_is_zero(arb_midref(x)) || prec == DECIMAL_PREC_EXACT)
        return _decball_set_arb_direct(res, x, ctx);

    fmpz_init(m);
    fmpz_init(t);
    _decmag_init(err, ctx);
    _decmag_init(eps, ctx);

    arf_get_fmpz_2exp(m, t, arb_midref(x));

    /* moderate sizes: exact conversion of the midpoint */
    r = fmpz_fits_si(t) ? _decfloat_set_round_fmpz_2exp_scaled(&res->mid, m, fmpz_get_si(t), prec, DECIMAL_CTX_RND(ctx), &info,
            DECIMAL_CTX_PRECISE_RADIUS(ctx) ? err : NULL, ctx) : 0;

    if (r == 1)
    {
        _decmag_set_mag(&res->rad, arb_radref(x), ctx);
        _decball_add_rounding_error(res, &info, err, prec, ctx);
        status = GR_SUCCESS;
        goto cleanup;
    }

    if (r == -1)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* huge sizes: an exactly representable midpoint */
    if (_decfloat_set_exact_decimal_fmpz_2exp(&res->mid, m, t, prec, ctx))
    {
        _decmag_set_mag(&res->rad, arb_radref(x), ctx);
        status = _decfloat_finalize(&res->mid, ctx);
        goto cleanup;
    }

    /* scale with a power of ten, at increasing precision until the
       scaling error is negligible compared to the radius */
    wp = _decimal_digits_to_bits(prec) + 30;
    wp_max = 64 * wp + 4096;

    for (;;)
    {
        status = _decball_set_arb_scaled(res, x, wp, eps, ctx);
        if (status != GR_SUCCESS)
            break;
        _decmag_mul_ui(eps, eps, 16384, ctx);
        if (wp >= wp_max || _decmag_cmp(eps, &res->rad, ctx) <= 0)
            break;
        wp *= 2;
    }

cleanup:
    fmpz_clear(m);
    fmpz_clear(t);
    _decmag_clear(err, ctx);
    _decmag_clear(eps, ctx);
    return status;
}

/* ------------------------------------------------------------------------- */
/*    Ball functions                                                         */
/* ------------------------------------------------------------------------- */

int
decball_via_arb(decball_t res, const decball_t x, int (*func)(arb_t, const arb_t, slong), gr_ctx_t ctx)
{
    arb_t a;
    slong wp;
    int status;

    wp = _decimal_digits_to_bits(DECIMAL_CTX_PREC(ctx));

    arb_init(a);
    status = decball_get_arb(a, x, wp, ctx);
    if (status == GR_SUCCESS)
        status = func(a, a, wp);
    if (status == GR_SUCCESS)
        status = decball_set_arb(res, a, ctx);
    arb_clear(a);
    return status;
}

int
decball_pow(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx)
{
    arb_t a, b;
    slong wp;
    int status;

    /* integer exponents */
    if (DECMAG_IS_ZERO(&y->rad) && _decfloat_is_int(&y->mid, ctx)
        && (DECFLOAT_IS_ZERO(&y->mid) || _decfloat_sci_exp_clamped(&y->mid, ctx) < 12))
    {
        fmpz_t n;
        fmpz_init(n);
        if (decfloat_get_fmpz(n, &y->mid, ctx) == GR_SUCCESS && fmpz_bits(n) < 40)
        {
            status = _decball_pow_fmpz_binexp(res, x, n, ctx);
            fmpz_clear(n);
            return status;
        }
        fmpz_clear(n);
    }

    wp = _decimal_digits_to_bits(DECIMAL_CTX_PREC(ctx));

    arb_init(a);
    arb_init(b);
    status = decball_get_arb(a, x, wp, ctx);
    status |= decball_get_arb(b, y, wp, ctx);

    if (status == GR_SUCCESS)
    {
        if (arb_is_negative(a))
            status = GR_DOMAIN;
        else if (!arb_is_nonnegative(a))
            status = GR_UNABLE;
        else
        {
            arb_pow(a, a, b, wp);
            if (arb_is_finite(a))
                status = decball_set_arb(res, a, ctx);
            else
                status = GR_UNABLE;
        }
    }

    arb_clear(a);
    arb_clear(b);
    return status;
}


/* ------------------------------------------------------------------------- */
/*    Small arguments and exact cases                                        */
/* ------------------------------------------------------------------------- */

/* Ziv's strategy cannot decide the rounding when the exact result lies
   extremely close to a representable number, which happens
   systematically in directed rounding modes for tiny arguments, e.g.
   sin(10^-1000000000) = 10^-1000000000 - (tiny). The helpers below
   handle f(x) = S + t where S is an exact decimal computed from the
   leading Taylor terms and t is a nonzero tail of known sign lying far
   below the rounding position. */

#define DECFLOAT_EXP_CLAMP (WORD_MAX / 16)

/* scientific exponent of a finite nonzero x, clamped to +/- EXP_CLAMP */
slong
_decfloat_sci_exp_clamped(const decfloat_t x, gr_ctx_t ctx)
{
    slong E;

    if (decfloat_get_sci_exp_si(&E, x, ctx))
        return FLINT_MAX(FLINT_MIN(E, DECFLOAT_EXP_CLAMP), -DECFLOAT_EXP_CLAMP);

    return (fmpz_sgn(&x->exp) > 0) ? DECFLOAT_EXP_CLAMP : -DECFLOAT_EXP_CLAMP;
}

/*
    Rounds S + t to prec digits where S is an exact, finite, nonzero
    decfloat and t is a nonzero real number with sign tail_sign and
    |t| < 10^tail_exp. Succeeds (returning 1) when tail_exp is at least
    prec + 4e digits below the leading digit of S, so that the tail
    only acts as a sticky bit in the rounding. Returns 0 if the
    precondition fails and -1 on an error (exponent limits).
*/
int
_decfloat_round_with_tail(decfloat_t res, const decfloat_t S, int tail_sign, slong tail_exp, slong prec, int rnd, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    slong ES, K;
    decfloat_t tiny;
    int status;

    if (DECFLOAT_IS_SPECIAL(S) || prec == DECIMAL_PREC_EXACT)
        return 0;

    ES = _decfloat_sci_exp_clamped(S, ctx);

    /* prec <= DECIMAL_PREC_MAX = WORD_MAX / 64, so this cannot overflow */
    if (tail_exp > ES - prec - 4 * e)
        return 0;

    /* a formal tiny term 10^(e (v - K - 3)) whose only role is to sit
       below the rounding horizon of _decfloat_add, which then treats it
       as a sticky bit ("eps") of the given sign */
    K = (prec + e - 1) / e + 4;
    decfloat_init(tiny, ctx);
    radix_integer_fit_limbs(&tiny->m, 1, DECIMAL_CTX_RADIX(ctx));
    tiny->m.d[0] = 1;
    tiny->m.size = 1;
    fmpz_sub_ui(&tiny->exp, &S->exp, K + 3);

    status = _decfloat_add(res, S, tiny, tail_sign < 0, prec, rnd, NULL, NULL, ctx);

    decfloat_clear(tiny, ctx);
    return (status == GR_SUCCESS) ? 1 : -1;
}

/* x^y where the result is a rational number that we can compute exactly
   (y = p/q with q a small power of 2 and 5, x^p a perfect q-th power);
   returns 1 if handled */
static int
_decfloat_pow_exact(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpq_t xq, yq, rq;
    fmpz_t n, d, t;
    ulong p, q;
    int r = 0, neg;

    fmpq_init(xq);
    fmpq_init(yq);
    fmpq_init(rq);
    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(t);

    /* y = p/q with |p|, q < 2^20 requires at most 7 digits and a decimal
       valuation of at least -6; x^p must be of moderate size */
    if (DECFLOAT_IS_SPECIAL(x) || DECFLOAT_IS_SPECIAL(y)
        || decfloat_digits(y, ctx) > 7 || _decfloat_val10_clamped(y, ctx) < -6
        || _decfloat_sci_exp_clamped(y, ctx) > 7
        || (double) (decfloat_digits(x, ctx) + FLINT_ABS(_decfloat_val10_clamped(x, ctx))) > 200000.0)
        goto cleanup;

    if (decfloat_get_fmpq(xq, x, ctx) != GR_SUCCESS || decfloat_get_fmpq(yq, y, ctx) != GR_SUCCESS)
        goto cleanup;

    if (!fmpz_fits_si(fmpq_numref(yq)) || fmpz_bits(fmpq_numref(yq)) > 20 || fmpz_bits(fmpq_denref(yq)) > 20)
        goto cleanup;

    p = fmpz_get_si(fmpq_numref(yq));
    neg = (fmpz_sgn(fmpq_numref(yq)) < 0);
    p = FLINT_ABS((slong) p);
    q = fmpz_get_ui(fmpq_denref(yq));

    if (p == 0)
        goto cleanup;

    /* avoid huge exact powers */
    if ((fmpz_bits(fmpq_numref(xq)) + fmpz_bits(fmpq_denref(xq))) * p > 200000)
        goto cleanup;

    fmpz_pow_ui(n, fmpq_numref(xq), p);
    fmpz_pow_ui(d, fmpq_denref(xq), p);

    if (q != 1)
    {
        if (fmpz_sgn(n) < 0)
            goto cleanup;
        fmpz_root(t, n, q);
        fmpz_pow_ui(t, t, q);
        if (!fmpz_equal(t, n))
            goto cleanup;
        fmpz_root(n, n, q);
        fmpz_root(t, d, q);
        fmpz_pow_ui(t, t, q);
        if (!fmpz_equal(t, d))
            goto cleanup;
        fmpz_root(d, d, q);
    }

    if (neg)
    {
        if (fmpz_is_zero(n))
            goto cleanup;
        fmpz_swap(n, d);
        if (fmpz_sgn(d) < 0)
        {
            fmpz_neg(d, d);
            fmpz_neg(n, n);
        }
    }

    fmpq_set_fmpz_frac(rq, n, d);
    r = (decfloat_set_round_fmpq(res, rq, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;

cleanup:
    fmpq_clear(xq);
    fmpq_clear(yq);
    fmpq_clear(rq);
    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(t);
    return r;
}

/* x^y = exp(y log x) with y log x tiny: 1 + (small tail); returns 1 if handled */
static int
_decfloat_pow_small(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
{
    decfloat_t one, h;
    slong Ex, Ey, Eh, L, tail_exp;
    int r = 0, tail_sign;

    decfloat_init(one, ctx);
    decfloat_init(h, ctx);
    decfloat_one(one, ctx);

    Ex = _decfloat_sci_exp_clamped(x, ctx);
    Ey = _decfloat_sci_exp_clamped(y, ctx);

    if (Ex == 0 || Ex == -1)
    {
        /* h = x - 1, exactly; only needed (and cheap to compute
           relative to the input) when x is close to 1 */
        if (_decfloat_add(h, x, one, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, ctx) != GR_SUCCESS)
            goto cleanup;

        if (DECFLOAT_IS_ZERO(h))
        {
            r = (decfloat_one(res, ctx) == GR_SUCCESS) ? 1 : -1;
            goto cleanup;
        }

        Eh = _decfloat_sci_exp_clamped(h, ctx);
    }
    else
    {
        /* |h| >= 9/10, and h has the sign of x - 1 */
        Eh = 0;
        GR_MUST_SUCCEED(decfloat_set_si(h, (Ex >= 0) ? 1 : -1, ctx));
    }

    /* |log x| < 10^L */
    if (Eh <= -1)
        L = Eh + 2;                    /* |log(1+h)| <= 2|h| for |h| <= 1/2 */
    else
    {
        slong a = FLINT_ABS(Ex) + 1;   /* |log x| < 3 (|E_x| + 1) */
        L = 1;
        while (a >= 10) { a /= 10; L++; }
        L += 1;
    }

    if (Ey > DECFLOAT_EXP_CLAMP / 2 || L > DECFLOAT_EXP_CLAMP / 2)
        goto cleanup;

    /* u = y log x, |u| < 10^(Ey + 1 + L); |e^u - 1| <= 2|u| for |u| <= 1 */
    tail_exp = Ey + L + 2;
    tail_sign = _decfloat_sgn(y, ctx) * _decfloat_sgn(h, ctx);

    r = _decfloat_round_with_tail(res, one, tail_sign, tail_exp, prec, rnd, ctx);

cleanup:
    decfloat_clear(one, ctx);
    decfloat_clear(h, ctx);
    return r;
}

/* x^n for an integer n, exactly when feasible; returns 1 if handled */
int
_decfloat_pow_int_exact(decfloat_t res, const decfloat_t x, const fmpz_t n, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz_t m, t;
    ulong an;
    int r = 0;

    fmpz_init(m);
    fmpz_init(t);

    /* powers of ten: (+/- 10^t)^n = (+/-1)^n 10^(t n) for any n */
    if (!DECFLOAT_IS_SPECIAL(x) && decfloat_digits(x, ctx) == 1
        && decfloat_get_digit_si(x, _decfloat_val10_clamped(x, ctx), ctx) == 1
        && decfloat_get_fmpz_10exp_fmpz(m, t, x, ctx) == GR_SUCCESS)
    {
        if (fmpz_sgn(m) < 0 && fmpz_is_even(n))
            fmpz_one(m);
        fmpz_mul(t, t, n);
        r = (decfloat_set_round_fmpz_10exp_fmpz(res, m, t, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
        fmpz_clear(m);
        fmpz_clear(t);
        return r;
    }

    if (!fmpz_fits_si(n))
    {
        fmpz_clear(m);
        fmpz_clear(t);
        return 0;
    }

    an = FLINT_ABS(fmpz_get_si(n));

    /* x = m 10^t with m not divisible by 10; the exact result has about
       bits(m) * |n| bits, which we cap at roughly 100000 digits */
    if (!DECFLOAT_IS_SPECIAL(x) && (double) decfloat_digits(x, ctx) * (double) an <= 100000.0
        && FLINT_ABS(_decfloat_val10_clamped(x, ctx)) < 1000000
        && decfloat_get_fmpz_10exp_fmpz(m, t, x, ctx) == GR_SUCCESS && fmpz_bits(t) <= 20)
    {
        fmpz_pow_ui(m, m, an);
        fmpz_mul_ui(t, t, an);

        if (fmpz_sgn(n) >= 0)
        {
            r = (decfloat_set_round_fmpz_10exp_fmpz(res, m, t, prec, rnd, ctx) == GR_SUCCESS) ? 1 : -1;
        }
        else
        {
            /* 1 / (m 10^t) = (1/m) 10^(-t): round 1/m, then shift exactly */
            fmpq_t q;
            fmpq_init(q);
            fmpz_one(fmpq_numref(q));
            fmpz_set(fmpq_denref(q), m);
            if (fmpz_sgn(fmpq_denref(q)) < 0)
            {
                fmpz_neg(fmpq_numref(q), fmpq_numref(q));
                fmpz_neg(fmpq_denref(q), fmpq_denref(q));
            }
            fmpz_neg(t, t);
            r = (decfloat_set_round_fmpq(res, q, prec, rnd | DECIMAL_RND_NOLIMITS, ctx) == GR_SUCCESS) ? 1 : -1;
            if (r == 1)
                r = (_decfloat_mul_10exp(res, res, t, prec, rnd, NULL, NULL, ctx) == GR_SUCCESS) ? 1 : -1;
            fmpq_clear(q);
        }
    }

    fmpz_clear(m);
    fmpz_clear(t);
    return r;
}

/* ------------------------------------------------------------------------- */
/*    Float functions (correctly rounded via Ziv's strategy)                 */
/* ------------------------------------------------------------------------- */

/* Given a ball Y = [mid +/- rad] (rad != 0), attempts to round to prec
   digits: succeeds if both endpoints round to the same value. */
int
_decfloat_round_ball(decfloat_t res, const decball_t Y, slong prec, int rnd, gr_ctx_t bctx, gr_ctx_t ctx)
{
    decfloat_t lo, hi, r, rl, rh;
    int status, ok = 0;

    if (DECMAG_IS_ZERO(&Y->rad))
        return decfloat_set_round(res, &Y->mid, prec, rnd, ctx) == GR_SUCCESS ? 1 : -1;

    if (DECMAG_IS_INF(&Y->rad) || DECFLOAT_IS_SPECIAL(&Y->mid))
        return 0;

    decfloat_init(lo, ctx);
    decfloat_init(hi, ctx);
    decfloat_init(r, ctx);
    decfloat_init(rl, ctx);
    decfloat_init(rh, ctx);

    status = _decmag_get_decfloat(r, &Y->rad, bctx);
    status |= _decfloat_add(lo, &Y->mid, r, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, bctx);
    status |= _decfloat_add(hi, &Y->mid, r, 0, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, bctx);
    status |= decfloat_set_round(rl, lo, prec, rnd, ctx);
    status |= decfloat_set_round(rh, hi, prec, rnd, ctx);

    if (status == GR_SUCCESS && decfloat_equal(rl, rh, ctx) == T_TRUE)
    {
        decfloat_swap(res, rl, ctx);
        ok = 1;
    }

    decfloat_clear(lo, ctx);
    decfloat_clear(hi, ctx);
    decfloat_clear(r, ctx);
    decfloat_clear(rl, ctx);
    decfloat_clear(rh, ctx);
    return ok;
}

/* Correctly rounded conversion of an arb ball (Ziv's strategy with the
   exact-scaling machinery of decball_set_arb). Returns GR_UNABLE if the ball
   is too wide to decide the rounding. If err is not NULL, it is set to a
   bound for the rounding error |x - res| which is tight (to about the
   radius precision) when x is exact, and zero when x is exact and
   exactly representable. */
int
_decfloat_set_round_arb_err(decfloat_t res, decmag_ptr err, const arb_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    gr_ctx_t bctx;
    decball_t Y;
    slong wp, wp_max;
    int status = GR_UNABLE, r, decided = 0;

    if (!arb_is_finite(x))
        return GR_UNABLE;

    if (arb_is_exact(x) && !arf_is_zero(arb_midref(x)))
    {
        fmpz_t m, t;
        fmpz_init(m);
        fmpz_init(t);
        arf_get_fmpz_2exp(m, t, arb_midref(x));
        r = _decfloat_set_exact_decimal_fmpz_2exp(res, m, t, prec, ctx);
        fmpz_clear(m);
        fmpz_clear(t);
        if (r)
        {
            if (err != NULL)
                _decmag_zero(err, ctx);
            return _decfloat_finalize(res, ctx);
        }
    }

    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec, DECIMAL_RND_DOWN, 0);
    decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);
    decball_init(Y, bctx);

    wp_max = 64 * prec + 1000;

    for (wp = prec + 10; wp <= wp_max; wp *= 2)
    {
        decimal_ctx_set_prec(bctx, wp);

        status = decball_set_arb(Y, x, bctx);
        if (status != GR_SUCCESS)
            break;

        r = _decfloat_round_ball(res, Y, prec, rnd, bctx, ctx);
        if (r == -1)
        {
            status = GR_UNABLE;
            break;
        }

        if (r == 1)
        {
            decided = 1;

            if (err == NULL)
            {
                status = _decfloat_finalize(res, ctx);
                break;
            }
            else
            {
                /* err = |mid - res| + rad; tight enough if rad is
                   negligible compared to |mid - res| */
                decfloat_t d;
                decmag_t dm, t;
                int tight;

                decfloat_init(d, bctx);
                _decmag_init(dm, bctx);
                _decmag_init(t, bctx);
                GR_MUST_SUCCEED(_decfloat_add(d, &Y->mid, res, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, NULL, NULL, bctx));
                _decmag_set_decfloat_lower(dm, d, bctx);
                _decmag_mul_ui(t, &Y->rad, 1000000000, bctx);
                tight = DECMAG_IS_ZERO(&Y->rad) || (_decmag_cmp(t, dm, bctx) <= 0);
                _decmag_set_decfloat(dm, d, bctx);
                _decmag_add(dm, dm, &Y->rad, bctx);
                /* convert to the radius precision of ctx */
                if (DECMAG_IS_ZERO(dm))
                    _decmag_zero(err, ctx);
                else if (DECMAG_IS_INF(dm))
                    _decmag_inf(err, ctx);
                else
                    _decmag_set_uiui_10exp_fmpz(err, 0, dm->m, 10, 0, &dm->exp, ctx);
                decfloat_clear(d, bctx);
                _decmag_clear(dm, bctx);
                _decmag_clear(t, bctx);

                if (tight || arb_is_exact(x) == 0 || wp * 2 > wp_max)
                {
                    status = _decfloat_finalize(res, ctx);
                    break;
                }
            }
        }

        status = GR_UNABLE;
    }

    if (decided && status == GR_UNABLE)
        status = _decfloat_finalize(res, ctx);

    decball_clear(Y, bctx);
    gr_ctx_clear(bctx);
    return status;
}

int
decfloat_via_arb(decfloat_t res, const decfloat_t x, int (*func)(arb_t, const arb_t, slong), gr_ctx_t ctx)
{
    slong prec = DECIMAL_CTX_PREC(ctx);
    int rnd = DECIMAL_CTX_RND(ctx);
    slong wp, wpbits;
    arb_t a;
    decball_t Y;
    gr_ctx_t bctx;
    int status = GR_UNABLE, r;

    if (DECFLOAT_IS_NAN(x))
        return decfloat_nan(res, ctx);

    if (prec == DECIMAL_PREC_EXACT)
        return GR_UNABLE;

    /* a ball context with the same radix for the intermediate conversions */
    _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec, DECIMAL_RND_DOWN, 0);
    decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);

    arb_init(a);
    decball_init(Y, bctx);

    for (wp = prec + 10; wp <= _decimal_ziv_wp_max(prec, decfloat_digits(x, ctx)); wp *= 2)
    {
        wpbits = _decimal_digits_to_bits(wp);

        decimal_ctx_set_prec(bctx, wp);

        status = decfloat_get_arb(a, x, wpbits, ctx);
        if (status == GR_SUCCESS)
            status = func(a, a, wpbits);
        if (status != GR_SUCCESS)
            break;

        if (!arb_is_finite(a))
        {
            if (arf_is_pos_inf(arb_midref(a)) && mag_is_finite(arb_radref(a)))
                status = decfloat_pos_inf(res, ctx);
            else if (arf_is_neg_inf(arb_midref(a)) && mag_is_finite(arb_radref(a)))
                status = decfloat_neg_inf(res, ctx);
            else
            {
                /* an infinite radius typically means an intermediate
                   overflow at this precision: try again with more */
                status = GR_UNABLE;
                continue;
            }
            break;
        }

        status = decball_set_arb(Y, a, bctx);
        if (status != GR_SUCCESS)
            break;

        r = _decfloat_round_ball(res, Y, prec, rnd, bctx, ctx);
        if (r == 1)
        {
            status = _decfloat_finalize(res, ctx);
            break;
        }
        if (r == -1)
        {
            status = GR_UNABLE;
            break;
        }

        status = GR_UNABLE;
    }

    arb_clear(a);
    decball_clear(Y, bctx);
    gr_ctx_clear(bctx);
    return status;
}

int
decfloat_pow(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
        return decfloat_nan(res, ctx);

    DECFLOAT_CHECK_OPERAND(x, ctx);
    DECFLOAT_CHECK_OPERAND(y, ctx);

    if (DECFLOAT_IS_ZERO(y))
        return decfloat_one(res, ctx);

    /* integer exponents below 10^1000; larger ones are treated as general
       real exponents, without converting them to an integer */
    if (_decfloat_is_int(y, ctx) && !DECFLOAT_IS_SPECIAL(y) && _decfloat_sci_exp_clamped(y, ctx) < 1000)
    {
        fmpz_t n;
        int status;
        fmpz_init(n);
        if (decfloat_get_fmpz(n, y, ctx) == GR_SUCCESS)
        {
            status = decfloat_pow_fmpz(res, x, n, ctx);
            fmpz_clear(n);
            return status;
        }
        fmpz_clear(n);
    }

    if (DECFLOAT_IS_ZERO(x))
    {
        if (_decfloat_sgn(y, ctx) > 0)
            return decfloat_zero(res, ctx);
        return GR_DOMAIN;
    }

    if (_decfloat_sgn(x, ctx) < 0)
        return GR_DOMAIN;

    if (DECFLOAT_IS_FINITE(x) && DECFLOAT_IS_FINITE(y))
    {
        int r = _decfloat_pow_exact(res, x, y, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
        if (r == 0 && !DECIMAL_CTX_IS_EXACT(ctx))
            r = _decfloat_pow_small(res, x, y, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
        if (r == 1) return GR_SUCCESS;
        if (r == -1) return GR_UNABLE;
    }

    if (DECIMAL_CTX_IS_EXACT(ctx))
        return GR_UNABLE;

    return _decfloat_pow_ziv(res, x, y, ctx);
}

/* x^y = exp(y log x): Ziv loop with two arguments; x > 0 finite, y finite */
int
_decfloat_pow_ziv(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    {
        slong prec = DECIMAL_CTX_PREC(ctx);
        int rnd = DECIMAL_CTX_RND(ctx);
        slong wp, wpbits;
        arb_t a, b;
        decball_t Z;
        gr_ctx_t bctx;
        int status = GR_UNABLE, r;

        _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec, DECIMAL_RND_DOWN, 0);
        decimal_ctx_set_rad_prec(bctx, DECMAG_MAX_PREC);

        arb_init(a);
        arb_init(b);
        decball_init(Z, bctx);

        for (wp = prec + 10; wp <= _decimal_ziv_wp_max(prec, FLINT_MAX(decfloat_digits(x, ctx), decfloat_digits(y, ctx))); wp *= 2)
        {
            wpbits = _decimal_digits_to_bits(wp);
            decimal_ctx_set_prec(bctx, wp);

            status = decfloat_get_arb(a, x, wpbits, ctx);
            status |= decfloat_get_arb(b, y, wpbits, ctx);
            if (status != GR_SUCCESS)
                break;

            arb_pow(a, a, b, wpbits);

            if (!arb_is_finite(a))
            {
                /* an infinite radius typically means an intermediate
                   overflow at this precision: try again with more */
                status = GR_UNABLE;
                continue;
            }

            status = decball_set_arb(Z, a, bctx);
            if (status != GR_SUCCESS)
                break;

            r = _decfloat_round_ball(res, Z, prec, rnd, bctx, ctx);
            if (r == 1)
            {
                status = _decfloat_finalize(res, ctx);
                break;
            }
            if (r == -1)
            {
                status = GR_UNABLE;
                break;
            }
            status = GR_UNABLE;
        }

        arb_clear(a);
        arb_clear(b);
        decball_clear(Z, bctx);
        gr_ctx_clear(bctx);
        return status;
    }
}

/* ------------------------------------------------------------------------- */
/*    Vectors                                                                */
/* ------------------------------------------------------------------------- */

int
decfloat_vec_dot(decfloat_t res, const decfloat_t initial, int subtract, decfloat_srcptr vec1, decfloat_srcptr vec2, slong len, gr_ctx_t ctx)
{
    return gr_generic_vec_dot(res, initial, subtract, vec1, vec2, len, ctx);
}

int
decfloat_vec_dot_rev(decfloat_t res, const decfloat_t initial, int subtract, decfloat_srcptr vec1, decfloat_srcptr vec2, slong len, gr_ctx_t ctx)
{
    return gr_generic_vec_dot_rev(res, initial, subtract, vec1, vec2, len, ctx);
}

POP_OPTIONS
