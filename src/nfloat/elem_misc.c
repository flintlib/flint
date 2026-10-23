/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "gr.h"
#include "elem.h"

/*
    The less common elementary functions, as compositions of the native
    functions and arithmetic operations in an inner nfloat context with
    m >= n limbs (one more limb than the output precision, or none when
    the function precision is a dozen bits below the output precision).

    Each composition comes with a bound K for its relative error in units
    of u = 2^(-FLINT_BITS m), counting
      4u for add, sub and mul (rounding error at most 2 ulp; 0 when exact
         by Sterbenz's lemma, as noted),
      2u for div, inv and sqrt (less than 1 ulp),
      3u for the native functions (at most 1 + 2^-10 ulp),
      4u for pi,
    and propagating errors to first order via relative condition numbers
    |x f'(x) / f(x)| (those of the functions used here being at most 1 for
    the arguments where they are applied). The result r is then converted
    to the output with an error of K + 2 ulp(r) (covering the second-order
    terms), giving valid bounds with directed rounding.

    On 64-bit machines, at one and two limbs, the reciprocal functions and
    sinc first try the fast paths of sin_cos.c and hyperbolic.c, and the
    inverse functions those of elem_arc_fast.c (the same compositions in
    fixed point without the inner context).

    Arguments whose evaluation over- or underflows in the inner context,
    and n = NFLOAT_MAX_LIMBS (no room for a guard limb), are evaluated
    with mp_real balls (_nfloat_elem_ball) at increasing precision.
    Exact special values (zeros, poles, multiples of 1/2 for the functions
    of pi x, |x| = 1 for the inverse functions) are handled before.
*/

typedef enum
{
    ELEM_EXP10, ELEM_LOG10,
    ELEM_COT, ELEM_SEC, ELEM_CSC, ELEM_SINC,
    ELEM_COT_PI, ELEM_SEC_PI, ELEM_CSC_PI, ELEM_SINC_PI,
    ELEM_COTH, ELEM_SECH, ELEM_CSCH,
    ELEM_ASIN, ELEM_ACOS, ELEM_ATAN, ELEM_ACOT, ELEM_ASEC, ELEM_ACSC,
    ELEM_ASINH, ELEM_ACOSH, ELEM_ATANH, ELEM_ACOTH, ELEM_ASECH, ELEM_ACSCH
}
elem_func_t;

/* error bounds K (see above), not including division by pi */
static const unsigned char _elem_err[] =
{
    0, 8,
    8, 5, 5, 5,
    8, 5, 5, 5,
    5, 5, 5,
    11, 11, 3, 3, 11, 11,
    23, 15, 9, 9, 17, 25
};

/* division by pi: 4u for pi, 2u for the division */
#define ELEM_ERR_DIV_PI 6

/* Ball evaluation (the slow path) ******************************************/

/* r = a / b, or 0 if b is not certified nonzero with a relative radius
   below 2^-32 (as mp_real_div requires) */
static int
_ball_div(mp_real_t r, const mp_real_t a, const mp_real_t b, slong n)
{
    if (b->size == 0 || mp_real_rel_radius_lt_2exp_si(b) > -32)
        return 0;
    mp_real_div(r, a, b, n);
    return 1;
}

/* r 2^shift = exp(t) for a ball t of magnitude below 2^(FLINT_BITS - 1)
   (t is destroyed): exp(t - k log 2) 2^k with k near t / log 2 for
   |t| >= 2^(FLINT_BITS - 8), the subtraction at enough precision for the
   cancellation (mp_real_exp_bits takes arguments below 2^(FLINT_BITS - 5)) */
static void
_ball_exp_2exp(mp_real_t r, slong * shift, mp_real_t t, slong prec)
{
    slong et = mp_real_abs_bound_lt_2exp_si(t);

    *shift = 0;

    if (et > FLINT_BITS - 8 && t->size != 0)
    {
        mp_real_t u, v;
        slong n2 = mp_real_prec_bits(prec + et + 32);
        double td;
        slong k;

        /* the top limb of the midpoint gives |t| within a relative 2^-31,
           so that |t - k log 2| < 2^(FLINT_BITS - 30) */
        td = ldexp((double) t->d[t->size - 1], FLINT_BITS * (t->exp - 1));
        k = (slong) (td * 1.4426950408889634);
        if (t->negative)
            k = -k;

        mp_real_init(u);
        mp_real_init(v);
        mp_real_const_log2(u, n2, 1);
        mp_real_set_si(v, k);
        mp_real_mul(u, u, v, n2);
        mp_real_sub(t, t, u, n2);
        mp_real_clear(u);
        mp_real_clear(v);
        *shift = k;
    }

    mp_real_exp_bits(r, t, prec);
}

/* f(x) or f(x, y) (op = 2 f + div_pi) at about prec bits, for the
   arguments that pass the special-value checks of the functions below;
   0 when prec does not suffice */
static int
_nfloat_elem_ball(mp_real_t r, slong * shift, const mp_real_t x, const mp_real_t y, int op, slong prec)
{
    elem_func_t f = (elem_func_t) (op >> 1);
    int div_pi = op & 1, ok = 1;
    slong n, n2, wp = prec + 8;
    mp_real_t a, b, c;

    n = mp_real_prec_bits(wp);

    mp_real_init(a);
    mp_real_init(b);
    mp_real_init(c);

    switch (f)
    {
        case ELEM_EXP10:
            /* t = x log(10) accurate to about 2^-wp absolutely */
            n2 = mp_real_prec_bits(wp + FLINT_MAX(mp_real_abs_bound_lt_2exp_si(x), 0) + 16);
            mp_real_const_log10(a, n2, 1);
            mp_real_mul(a, a, x, n2);
            _ball_exp_2exp(r, shift, a, wp);
            break;

        case ELEM_LOG10:
            if ((ok = mp_real_log_bits(a, x, wp)))
            {
                mp_real_const_log10(b, n, 1);
                ok = _ball_div(r, a, b, n);
            }
            break;

        case ELEM_COT: case ELEM_COT_PI:
            /* 1 / tan, which keeps its relative accuracy near the zeros
               and poles */
            ok = (f == ELEM_COT) ? mp_real_tan_bits(a, x, wp) : mp_real_tan_pi_bits(a, x, wp);
            if (ok)
            {
                mp_real_set_ui(b, 1);
                ok = _ball_div(r, b, a, n);
            }
            break;

        case ELEM_SEC: case ELEM_CSC: case ELEM_SINC:
        case ELEM_SEC_PI: case ELEM_CSC_PI: case ELEM_SINC_PI:
            if (f >= ELEM_COT_PI)
                mp_real_sin_cos_pi_bits(a, b, x, wp);
            else
                mp_real_sin_cos_bits(a, b, x, wp);
            if (f == ELEM_SINC || f == ELEM_SINC_PI)
            {
                ok = _ball_div(r, a, x, n);
                if (f == ELEM_SINC_PI)
                    div_pi = 1;
            }
            else
            {
                mp_real_set_ui(c, 1);
                ok = _ball_div(r, c, (f == ELEM_SEC || f == ELEM_SEC_PI) ? b : a, n);
            }
            break;

        case ELEM_COTH:
            mp_real_tanh_bits(a, x, wp);
            mp_real_set_ui(b, 1);
            ok = _ball_div(r, b, a, n);
            break;

        case ELEM_SECH:
        case ELEM_CSCH:
            if (mp_real_abs_bound_lt_2exp_si(x) > FLINT_BITS - 6)
            {
                /* sech x = 2 exp(-|x|) (1 - delta), |csch x| =
                   2 exp(-|x|) (1 + delta'), 0 < delta, delta' <
                   4 exp(-2|x|) < 2^(-2^(FLINT_BITS - 7)), far below
                   2^-(wp + 64) */
                mp_real_set(a, x);
                a->negative = 1;
                _ball_exp_2exp(r, shift, a, wp);
                *shift += 1;
                mp_real_add_rel_error_2exp_si(r, -wp - 64);
                if (f == ELEM_CSCH)
                    r->negative = x->negative;
            }
            else
            {
                if (f == ELEM_SECH)
                    mp_real_sinh_cosh_bits(NULL, a, x, wp);
                else
                    mp_real_sinh_cosh_bits(a, NULL, x, wp);
                mp_real_set_ui(b, 1);
                ok = _ball_div(r, b, a, n);
            }
            break;

        case ELEM_ASIN: ok = mp_real_asin_bits(r, x, wp); break;
        case ELEM_ACOS: ok = mp_real_acos_bits(r, x, wp); break;
        case ELEM_ATAN: mp_real_atan_bits(r, x, wp); break;
        case ELEM_ASINH: mp_real_asinh_bits(r, x, wp); break;
        case ELEM_ACOSH: ok = mp_real_acosh_bits(r, x, wp); break;
        case ELEM_ATANH: ok = mp_real_atanh_bits(r, x, wp); break;

        case ELEM_ACOT: case ELEM_ASEC: case ELEM_ACSC:
        case ELEM_ACOTH: case ELEM_ASECH: case ELEM_ACSCH:
            /* the inverse functions at 1/x (a ball, so that next to the
               boundaries |x| = 1 the domain checks fail until the
               precision suffices) */
            mp_real_set_ui(b, 1);
            mp_real_div(a, b, x, n);
            switch (f)
            {
                case ELEM_ACOT: mp_real_atan_bits(r, a, wp); break;
                case ELEM_ASEC: ok = mp_real_acos_bits(r, a, wp); break;
                case ELEM_ACSC: ok = mp_real_asin_bits(r, a, wp); break;
                case ELEM_ACOTH: ok = mp_real_atanh_bits(r, a, wp); break;
                case ELEM_ASECH: ok = mp_real_acosh_bits(r, a, wp); break;
                default: mp_real_asinh_bits(r, a, wp); break;
            }
            break;
    }

    if (ok && div_pi)
    {
        mp_real_const_pi4(a, n, 1);
        mp_real_mul_2exp_si(a, a, 2);
        mp_real_div(r, r, a, n);
    }

    mp_real_clear(a);
    mp_real_clear(b);
    mp_real_clear(c);
    return ok;
}

/* y (in ictx, m >= n limbs) = x (in ctx, n limbs) exactly; x normal */
static void
_elem_get(nfloat_ptr y, nfloat_srcptr x, gr_ctx_t ictx, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx), m = NFLOAT_CTX_NLIMBS(ictx);

    NFLOAT_EXP(y) = NFLOAT_EXP(x);
    NFLOAT_SGNBIT(y) = NFLOAT_SGNBIT(x);
    _nfloat_zero_limbs(NFLOAT_D(y), m - n);
    _nfloat_copy_limbs(NFLOAT_D(y) + m - n, NFLOAT_D(x), n);
}

/* (the operations leave garbage on over- or underflow) */
#define CHK(op) do { if ((op) != GR_SUCCESS) return GR_UNABLE; } while (0)

/* evaluates f(x) (or f(x, y)), optionally divided by pi, in an inner
   context; GR_UNABLE if something over- or underflows there */
static int
_nfloat_elem_inner(nfloat_ptr res, elem_func_t f, int div_pi, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    ulong X[NFLOAT_MAX_ALLOC], A[NFLOAT_MAX_ALLOC], B[NFLOAT_MAX_ALLOC];
    ulong S[NFLOAT_MAX_ALLOC], R[NFLOAT_MAX_ALLOC], O[NFLOAT_MAX_ALLOC];
    ulong t[NFLOAT_MAX_LIMBS + 1];
    gr_ctx_t ictx;
    slong m, e = 0;
    ulong K;
    int sgn = NFLOAT_SGNBIT(x);

    K = _elem_err[f] + (div_pi ? ELEM_ERR_DIV_PI : 0);

    if (f == ELEM_EXP10)
    {
        /* t = x log(10) has relative error 7u (3u for log(10), 4u for the
           multiplication), i.e. absolute error 7u 2^e where |t| < 2^e;
           exp(t) adds 3u */
        e = FLINT_MAX(NFLOAT_EXP(x) + 2, 0);
        K = 7 * (UWORD(1) << e) + 3;
    }

    /* m limbs: the function precision plus 10 bits beyond the error */
    m = (NFLOAT_CTX_FUNC_PREC(ctx) + 10 + FLINT_BIT_COUNT(K + 2) + FLINT_BITS - 1) / FLINT_BITS;
    m = FLINT_MAX(m, NFLOAT_CTX_NLIMBS(ctx));
    if (m > NFLOAT_MAX_LIMBS)
        return GR_UNABLE;

    nfloat_ctx_init(ictx, FLINT_BITS * m, 0);
    nfloat_one(O, ictx);

    /* X = |x| */
    _elem_get(X, x, ictx, ctx);
    NFLOAT_SGNBIT(X) = 0;

    switch (f)
    {
        case ELEM_EXP10:
            NFLOAT_SGNBIT(X) = sgn;
            CHK(nfloat_set_ui(A, 10, ictx));
            CHK(nfloat_log(A, A, ictx));
            CHK(nfloat_mul(A, A, X, ictx));
            CHK(nfloat_exp(R, A, ictx));
            break;

        case ELEM_LOG10:
            CHK(nfloat_log(A, X, ictx));
            CHK(nfloat_set_ui(B, 10, ictx));
            CHK(nfloat_log(B, B, ictx));
            CHK(nfloat_div(R, A, B, ictx));
            break;

        case ELEM_COT: case ELEM_SEC: case ELEM_CSC: case ELEM_SINC:
        case ELEM_COT_PI: case ELEM_SEC_PI: case ELEM_CSC_PI: case ELEM_SINC_PI:
            NFLOAT_SGNBIT(X) = sgn;
            if (f >= ELEM_COT_PI)
                CHK(nfloat_sin_cos_pi(S, A, X, ictx));
            else
                CHK(nfloat_sin_cos(S, A, X, ictx));

            if (f == ELEM_COT || f == ELEM_COT_PI)
                CHK(nfloat_div(R, A, S, ictx));
            else if (f == ELEM_SEC || f == ELEM_SEC_PI)
                CHK(nfloat_inv(R, A, ictx));
            else if (f == ELEM_CSC || f == ELEM_CSC_PI)
                CHK(nfloat_inv(R, S, ictx));
            else
            {
                CHK(nfloat_div(R, S, X, ictx));
                if (f == ELEM_SINC_PI)
                    div_pi = 1;
            }
            break;

        case ELEM_COTH:
            NFLOAT_SGNBIT(X) = sgn;
            CHK(nfloat_tanh(A, X, ictx));
            CHK(nfloat_inv(R, A, ictx));
            break;

        case ELEM_SECH:
            CHK(nfloat_cosh(A, X, ictx));
            CHK(nfloat_inv(R, A, ictx));
            break;

        case ELEM_CSCH:
            NFLOAT_SGNBIT(X) = sgn;
            CHK(nfloat_sinh(A, X, ictx));
            CHK(nfloat_inv(R, A, ictx));
            break;

        case ELEM_ASIN:
        case ELEM_ACOS:
        case ELEM_ASEC:
        case ELEM_ACSC:
            /* S = sqrt((1 - |x|)(1 + |x|)) or sqrt((|x| - 1)(|x| + 1)),
               8u (the difference is exact for 1/2 <= |x| <= 2) */
            if (f == ELEM_ASIN || f == ELEM_ACOS)
                CHK(nfloat_sub(A, O, X, ictx));
            else
                CHK(nfloat_sub(A, X, O, ictx));
            CHK(nfloat_add(B, O, X, ictx));
            CHK(nfloat_mul(A, A, B, ictx));
            CHK(nfloat_sqrt(S, A, ictx));

            /* asin x = atan2(x, S), acos x = atan2(S, x),
               asec x = atan2(S, sgn(x)), acsc x = atan2(sgn(x), S) */
            NFLOAT_SGNBIT(X) = sgn;
            if (f == ELEM_ASEC || f == ELEM_ACSC)
                CHK(nfloat_set_si(X, sgn ? -1 : 1, ictx));
            if (f == ELEM_ASIN || f == ELEM_ACSC)
                CHK(nfloat_atan2(R, X, S, ictx));
            else
                CHK(nfloat_atan2(R, S, X, ictx));
            break;

        case ELEM_ATAN:
            NFLOAT_SGNBIT(X) = sgn;
            CHK(nfloat_atan(R, X, ictx));
            break;

        case ELEM_ACOT:
            /* acot x = atan2(sgn(x), |x|) */
            CHK(nfloat_set_si(A, sgn ? -1 : 1, ictx));
            CHK(nfloat_atan2(R, A, X, ictx));
            break;

        case ELEM_ACSCH:
            /* asinh(1/|x|), 2u + 23u */
            CHK(nfloat_inv(X, X, ictx));
            FLINT_FALLTHROUGH;

        case ELEM_ASINH:
            /* log1p(y + y^2 / (1 + sqrt(1 + y^2))), y = |x| */
            CHK(nfloat_sqr(A, X, ictx));                  /* 4u */
            CHK(nfloat_add(B, A, O, ictx));            /* 8u */
            CHK(nfloat_sqrt(B, B, ictx));                 /* 6u */
            CHK(nfloat_add(B, B, O, ictx));            /* 10u */
            CHK(nfloat_div(A, A, B, ictx));               /* 16u */
            CHK(nfloat_add(A, A, X, ictx));               /* 20u */
            CHK(nfloat_log1p(R, A, ictx));                /* 23u */
            NFLOAT_SGNBIT(R) = sgn;
            break;

        case ELEM_ACOSH:
            /* log1p(d + sqrt(d (x + 1))), d = x - 1 */
            CHK(nfloat_sub(A, X, O, ictx));            /* 4u, exact for x <= 2 */
            CHK(nfloat_add(B, X, O, ictx));            /* 4u */
            CHK(nfloat_mul(B, A, B, ictx));               /* 12u */
            CHK(nfloat_sqrt(B, B, ictx));                 /* 8u */
            CHK(nfloat_add(A, A, B, ictx));               /* 12u */
            CHK(nfloat_log1p(R, A, ictx));                /* 15u */
            break;

        case ELEM_ATANH:
        case ELEM_ACOTH:
            /* log1p(2|x| / (1 - |x|)) / 2, log1p(2 / (|x| - 1)) / 2 */
            if (f == ELEM_ATANH)
            {
                CHK(nfloat_sub(A, O, X, ictx));           /* 4u, exact for |x| >= 1/2 */
                CHK(nfloat_div(A, X, A, ictx));           /* 6u */
            }
            else
            {
                CHK(nfloat_sub(A, X, O, ictx));        /* 4u, exact for |x| <= 2 */
                CHK(nfloat_inv(A, A, ictx));              /* 6u */
            }
            NFLOAT_EXP(A) += 1;
            CHK(nfloat_log1p(R, A, ictx));                /* 9u */
            NFLOAT_EXP(R) -= 1;
            NFLOAT_SGNBIT(R) = sgn;
            break;

        case ELEM_ASECH:
            /* log1p((a + sqrt(a (1 + x))) / x), a = 1 - x */
            CHK(nfloat_sub(A, O, X, ictx));               /* 4u, exact for x >= 1/2 */
            CHK(nfloat_add(B, X, O, ictx));            /* 4u */
            CHK(nfloat_mul(B, A, B, ictx));               /* 12u */
            CHK(nfloat_sqrt(B, B, ictx));                 /* 8u */
            CHK(nfloat_add(A, A, B, ictx));               /* 12u */
            CHK(nfloat_div(A, A, X, ictx));               /* 14u */
            CHK(nfloat_log1p(R, A, ictx));                /* 17u */
            break;
    }

    if (div_pi)
    {
        CHK(nfloat_pi(A, ictx));
        CHK(nfloat_div(R, R, A, ictx));
    }

    if (NFLOAT_IS_SPECIAL(R))
        return GR_UNABLE;

    _nfloat_copy_limbs(t, NFLOAT_D(R), m);
    return _nfloat_set_mpn_err(res, t, m, NFLOAT_EXP(R) - FLINT_BITS * m, K + 2, K + 2, NFLOAT_SGNBIT(R), ctx);
}

static int
_nfloat_elem(nfloat_ptr res, elem_func_t f, int div_pi, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    int status;

#if FLINT_BITS == 64
    /* the inverse functions at one and two limbs (elem_arc_fast.c) */
    if (f >= ELEM_ASIN && f <= ELEM_ACSCH && NFLOAT_CTX_NLIMBS(ctx) <= 2)
    {
        static const unsigned char arc[] = {
            NFLOAT_ARC_ASIN, NFLOAT_ARC_ACOS, NFLOAT_ARC_ATAN, NFLOAT_ARC_ACOT,
            NFLOAT_ARC_ASEC, NFLOAT_ARC_ACSC, NFLOAT_ARC_ASINH, NFLOAT_ARC_ACOSH,
            NFLOAT_ARC_ATANH, NFLOAT_ARC_ACOTH, NFLOAT_ARC_ASECH, NFLOAT_ARC_ACSCH };

        status = _nfloat_arc_fast(res, x, arc[f - ELEM_ASIN], div_pi, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif

    status = _nfloat_elem_inner(res, f, div_pi, x, y, ctx);

    if (status == GR_UNABLE)
        status = _nfloat_eval_ball(res, _nfloat_elem_ball, x, y, 2 * f + div_pi,
            (f >= ELEM_COT && f <= ELEM_SINC) ? _nfloat_trig_ball_min_prec(x) : 0, ctx);

    return status;
}

/* the value c pi/4 (or c/4 for the functions divided by pi), sign sgnbit */
static int
_nfloat_pi4_mul(nfloat_ptr res, ulong c, int div_pi, int sgnbit, gr_ctx_t ctx)
{
    if (div_pi)
    {
        int status = nfloat_set_ui(res, c, ctx);
        NFLOAT_EXP(res) -= 2;
        NFLOAT_SGNBIT(res) = sgnbit;
        return status;
    }
    else
    {
        ulong t[NFLOAT_ELEM_MAX_LIMBS + 2];
        slong w = _nfloat_elem_limbs(ctx, 0);
        t[w] = mpn_mul_1(t, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w), w, c);
        return _nfloat_set_mpn_err(res, t, w + 1, -FLINT_BITS * w, c, c, sgnbit, ctx);
    }
}

/* |x| compared with 1: -1, 0, 1 */
static int
_nfloat_cmpabs_one(nfloat_srcptr x, slong n)
{
    if (NFLOAT_EXP(x) <= 0)
        return -1;
    if (NFLOAT_EXP(x) >= 2)
        return 1;
    if (NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(NFLOAT_D(x), n - 1))
        return 0;
    return 1;
}

/* exponential and logarithm in base 10 ************************************/

int
nfloat_exp10(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_exp(res, x, ctx);

    /* 10^x over- or underflows for |x| >= 2^(FLINT_BITS - 4):
       2^(FLINT_BITS - 4) log2(10) > 2^(FLINT_BITS - 2) */
    if (NFLOAT_EXP(x) > FLINT_BITS - 4)
        return NFLOAT_SGNBIT(x) ? _nfloat_underflow(res, 0, ctx) : _nfloat_overflow(res, 0, ctx);

    /* (bounds the error in _nfloat_elem_inner) */
    if (NFLOAT_EXP(x) > FLINT_BITS - 7)
        return _nfloat_eval_ball(res, _nfloat_elem_ball, x, NULL, 2 * ELEM_EXP10, 0, ctx);

    /* 10^x = 1 + delta, |delta| < 2.4 |x| */
    if (NFLOAT_EXP(x) < -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx) - 2)
        return _nfloat_one_plus_tiny(res, 0, NFLOAT_SGNBIT(x), -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx), ctx);

    return _nfloat_elem(res, ELEM_EXP10, 0, x, NULL, ctx);
}

int
nfloat_log10(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_SGNBIT(x))
        return nfloat_log(res, x, ctx);

    if (_nfloat_cmpabs_one(x, NFLOAT_CTX_NLIMBS(ctx)) == 0)
        return nfloat_zero(res, ctx);

    return _nfloat_elem(res, ELEM_LOG10, 0, x, NULL, ctx);
}

/* reciprocal trigonometric and hyperbolic functions *************************/

int
nfloat_cot(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_nan(res, ctx);
#if FLINT_BITS == 64
    {
        int status = _nfloat_trig_recip_fast(res, x, 0, NFLOAT_TRIG_COT, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif
    return _nfloat_elem(res, ELEM_COT, 0, x, NULL, ctx);
}

int
nfloat_sec(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return NFLOAT_IS_ZERO(x) ? nfloat_one(res, ctx) : nfloat_nan(res, ctx);
    /* 1 + delta, 0 < delta < x^2 */
    if (2 * NFLOAT_EXP(x) <= -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx))
        return _nfloat_one_plus_tiny(res, 0, 0, -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx), ctx);
#if FLINT_BITS == 64
    {
        int status = _nfloat_trig_recip_fast(res, x, 0, NFLOAT_TRIG_SEC, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif
    return _nfloat_elem(res, ELEM_SEC, 0, x, NULL, ctx);
}

int
nfloat_csc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_nan(res, ctx);
#if FLINT_BITS == 64
    {
        int status = _nfloat_trig_recip_fast(res, x, 0, NFLOAT_TRIG_CSC, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif
    return _nfloat_elem(res, ELEM_CSC, 0, x, NULL, ctx);
}

int
nfloat_sinc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return NFLOAT_IS_ZERO(x) ? nfloat_one(res, ctx) : nfloat_nan(res, ctx);
    /* 1 - delta, delta < x^2/6 */
    if (2 * NFLOAT_EXP(x) <= -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx))
        return _nfloat_one_plus_tiny(res, 0, 1, -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx), ctx);
#if FLINT_BITS == 64
    {
        int status = _nfloat_trig_recip_fast(res, x, 0, NFLOAT_TRIG_SINC, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif
    return _nfloat_elem(res, ELEM_SINC, 0, x, NULL, ctx);
}

int
nfloat_coth(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_nan(res, ctx);
#if FLINT_BITS == 64
    {
        int status = _nfloat_tanh_coth_fast(res, x, 1, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif

    return _nfloat_elem(res, ELEM_COTH, 0, x, NULL, ctx);
}

int
nfloat_sech(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return NFLOAT_IS_ZERO(x) ? nfloat_one(res, ctx) : nfloat_nan(res, ctx);
    /* sech x < 2 exp(-|x|) underflows for |x| >= 2^(FLINT_BITS - 3) */
    if (NFLOAT_EXP(x) > FLINT_BITS - 3)
        return _nfloat_underflow(res, 0, ctx);
    /* 1 - delta, 0 < delta < x^2 / 2 */
    if (2 * NFLOAT_EXP(x) <= -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx))
        return _nfloat_one_plus_tiny(res, 0, 1, -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx), ctx);
#if FLINT_BITS == 64
    {
        int status = _nfloat_tanh_coth_fast(res, x, 2, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif
    /* (cosh overflows in the inner context) */
    if (NFLOAT_EXP(x) > FLINT_BITS - 6)
        return _nfloat_eval_ball(res, _nfloat_elem_ball, x, NULL, 2 * ELEM_SECH, 0, ctx);
    return _nfloat_elem(res, ELEM_SECH, 0, x, NULL, ctx);
}

int
nfloat_csch(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_nan(res, ctx);
    if (NFLOAT_EXP(x) > FLINT_BITS - 3)
        return _nfloat_underflow(res, NFLOAT_SGNBIT(x), ctx);
#if FLINT_BITS == 64
    {
        int status = _nfloat_tanh_coth_fast(res, x, 3, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif
    if (NFLOAT_EXP(x) > FLINT_BITS - 6)
        return _nfloat_eval_ball(res, _nfloat_elem_ball, x, NULL, 2 * ELEM_CSCH, 0, ctx);
    return _nfloat_elem(res, ELEM_CSCH, 0, x, NULL, ctx);
}

/* functions of pi x */

static int
_nfloat_trig_pi_recip(nfloat_ptr res, nfloat_srcptr x, elem_func_t f, gr_ctx_t ctx)
{
    ulong T[NFLOAT_MAX_LIMBS + 1], xbuf[NFLOAT_MAX_LIMBS + 1];
    mp_real_t X;
    slong tl, te;
    int k;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x) && (f == ELEM_SEC_PI || f == ELEM_SINC_PI))
            return nfloat_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

#if FLINT_BITS == 64
    {
        int status = _nfloat_trig_recip_fast(res, x, 1,
            (f == ELEM_COT_PI) ? NFLOAT_TRIG_COT : (f == ELEM_SEC_PI) ? NFLOAT_TRIG_SEC :
            (f == ELEM_CSC_PI) ? NFLOAT_TRIG_CSC : NFLOAT_TRIG_SINC, ctx);
        if (status != NFLOAT_FAST_DECLINE)
            return status;
    }
#endif

    /* the exact reduction |x| = a/2 + s t, t in [0, 1/4]
       (mp_real_sin_cos_pi_bits); t = 0 for the multiples of 1/2 */
    _nfloat_mp_real_view(X, xbuf, x, ctx);
    k = _mp_real_trig_pi_reduce(T, &tl, &te, X) & 3;

    if (tl == 0)
    {
        /* x = k/2 mod 2: sin(pi x) = 0, 1, 0, -1; cos(pi x) = 1, 0, -1, 0 */
        if (NFLOAT_SGNBIT(x))
            k = (4 - k) & 3;

        switch (f)
        {
            case ELEM_COT_PI:
                return (k & 1) ? nfloat_zero(res, ctx) : nfloat_nan(res, ctx);
            case ELEM_SEC_PI:
                return (k & 1) ? nfloat_nan(res, ctx) : nfloat_set_si(res, (k == 0) ? 1 : -1, ctx);
            case ELEM_CSC_PI:
                return (k & 1) ? nfloat_set_si(res, (k == 1) ? 1 : -1, ctx) : nfloat_nan(res, ctx);
            default:
                /* sinc_pi: sin(pi x) / (pi x), x != 0 */
                if (!(k & 1))
                    return nfloat_zero(res, ctx);
                break;
        }
    }

    /* sinc_pi of tiny x: 1 - delta, delta < (pi x)^2 / 6 */
    if (f == ELEM_SINC_PI && 2 * NFLOAT_EXP(x) + 4 <= -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx))
        return _nfloat_one_plus_tiny(res, 0, 1, -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx), ctx);

    /* (the compositions reduce x exactly again) */
    return _nfloat_elem(res, f, 0, x, NULL, ctx);
}

int nfloat_cot_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_trig_pi_recip(res, x, ELEM_COT_PI, ctx); }
int nfloat_sec_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_trig_pi_recip(res, x, ELEM_SEC_PI, ctx); }
int nfloat_csc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_trig_pi_recip(res, x, ELEM_CSC_PI, ctx); }
int nfloat_sinc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_trig_pi_recip(res, x, ELEM_SINC_PI, ctx); }

/* inverse trigonometric functions ******************************************/

/* asin (which = 0) or acos (which = 1) of x, divided by pi if div_pi */
static int
_nfloat_asin_acos(nfloat_ptr res, nfloat_srcptr x, int which, int div_pi, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    int c;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return (which == 0) ? nfloat_zero(res, ctx) : _nfloat_pi4_mul(res, 2, div_pi, 0, ctx);
        return nfloat_nan(res, ctx);
    }

    c = _nfloat_cmpabs_one(x, n);

    if (c > 0)
        return nfloat_nan(res, ctx);

    if (c == 0)
    {
        /* asin(+-1) = +-pi/2, acos(1) = 0, acos(-1) = pi */
        if (which == 0)
            return _nfloat_pi4_mul(res, 2, div_pi, NFLOAT_SGNBIT(x), ctx);
        if (NFLOAT_SGNBIT(x))
            return _nfloat_pi4_mul(res, 4, div_pi, 0, ctx);
        return nfloat_zero(res, ctx);
    }

    /* asin x = x (1 + delta), delta < x^2 / 6 */
    if (which == 0 && !div_pi && 2 * NFLOAT_EXP(x) <= -FLINT_BITS * n)
        return _nfloat_set_tiny_rel_err(res, x, 1, ctx);

    return _nfloat_elem(res, (which == 0) ? ELEM_ASIN : ELEM_ACOS, div_pi, x, NULL, ctx);
}

int nfloat_asin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_asin_acos(res, x, 0, 0, ctx); }
int nfloat_acos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_asin_acos(res, x, 1, 0, ctx); }
int nfloat_asin_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_asin_acos(res, x, 0, 1, ctx); }
int nfloat_acos_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_asin_acos(res, x, 1, 1, ctx); }

int
nfloat_atan_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_INF(x))
            return _nfloat_pi4_mul(res, 2, 1, NFLOAT_IS_NEG_INF(x), ctx);
        return nfloat_nan(res, ctx);
    }

    if (_nfloat_cmpabs_one(x, NFLOAT_CTX_NLIMBS(ctx)) == 0)
        return _nfloat_pi4_mul(res, 1, 1, NFLOAT_SGNBIT(x), ctx);

    return _nfloat_elem(res, ELEM_ATAN, 1, x, NULL, ctx);
}

/* the inverse functions of 1/x */
static int
_nfloat_arc_recip(nfloat_ptr res, nfloat_srcptr x, elem_func_t f, int div_pi, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    int c;

    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_nan(res, ctx);

    c = _nfloat_cmpabs_one(x, n);

    switch (f)
    {
        case ELEM_ACOT:
            /* acot(+-1) = +-pi/4 */
            if (c == 0)
                return _nfloat_pi4_mul(res, 1, div_pi, NFLOAT_SGNBIT(x), ctx);
            break;
        case ELEM_ASEC:
            if (c < 0)
                return nfloat_nan(res, ctx);
            if (c == 0)
                return NFLOAT_SGNBIT(x) ? _nfloat_pi4_mul(res, 4, div_pi, 0, ctx) : nfloat_zero(res, ctx);
            break;
        case ELEM_ACSC:
            if (c < 0)
                return nfloat_nan(res, ctx);
            if (c == 0)
                return _nfloat_pi4_mul(res, 2, div_pi, NFLOAT_SGNBIT(x), ctx);
            break;
        case ELEM_ACOTH:
            if (c <= 0)
                return nfloat_nan(res, ctx);
            break;
        case ELEM_ASECH:
            if (NFLOAT_SGNBIT(x) || c > 0)
                return nfloat_nan(res, ctx);
            if (c == 0)
                return nfloat_zero(res, ctx);
            break;
        default:
            break;
    }

    return _nfloat_elem(res, f, div_pi, x, NULL, ctx);
}

int nfloat_acot(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ACOT, 0, ctx); }
int nfloat_asec(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ASEC, 0, ctx); }
int nfloat_acsc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ACSC, 0, ctx); }
int nfloat_acot_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ACOT, 1, ctx); }
int nfloat_asec_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ASEC, 1, ctx); }
int nfloat_acsc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ACSC, 1, ctx); }
int nfloat_acoth(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ACOTH, 0, ctx); }
int nfloat_asech(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ASECH, 0, ctx); }
int nfloat_acsch(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx) { return _nfloat_arc_recip(res, x, ELEM_ACSCH, 0, ctx); }

/* inverse hyperbolic functions *********************************************/

int
nfloat_asinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_set(res, x, ctx) | (NFLOAT_IS_NAN(x) ? nfloat_nan(res, ctx) : GR_SUCCESS);

    /* x (1 - delta), delta < x^2 / 6 */
    if (2 * NFLOAT_EXP(x) <= -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx))
        return _nfloat_set_tiny_rel_err(res, x, -1, ctx);

    return _nfloat_elem(res, ELEM_ASINH, 0, x, NULL, ctx);
}

int
nfloat_acosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    int c;

    if (NFLOAT_IS_SPECIAL(x))
        return NFLOAT_IS_POS_INF(x) ? nfloat_pos_inf(res, ctx) : nfloat_nan(res, ctx);

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    c = _nfloat_cmpabs_one(x, NFLOAT_CTX_NLIMBS(ctx));
    if (c < 0)
        return nfloat_nan(res, ctx);
    if (c == 0)
        return nfloat_zero(res, ctx);

    return _nfloat_elem(res, ELEM_ACOSH, 0, x, NULL, ctx);
}

int
nfloat_atanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return NFLOAT_IS_ZERO(x) ? nfloat_zero(res, ctx) : nfloat_nan(res, ctx);

    if (_nfloat_cmpabs_one(x, NFLOAT_CTX_NLIMBS(ctx)) >= 0)
        return nfloat_nan(res, ctx);

    /* x (1 + delta), delta < x^2 / 3 / (1 - x^2) */
    if (2 * NFLOAT_EXP(x) <= -FLINT_BITS * NFLOAT_CTX_NLIMBS(ctx))
        return _nfloat_set_tiny_rel_err(res, x, 1, ctx);

    return _nfloat_elem(res, ELEM_ATANH, 0, x, NULL, ctx);
}

/*
    hypot of normal x, y with |x| >= |y| and d = EXP(x) - EXP(y) <=
    FLINT_BITS / 2 n + 2, exactly: with the n-limb mantissas X and Y,
    x^2 + y^2 = S 2^E for the integer S = X^2 2^(2d) + Y^2 (at most
    3n + 2 limbs) and the even exponent E = 2 EXP(x) - 2 FLINT_BITS n - 2d,
    so sqrt(x^2 + y^2) = sqrt(S) 2^(E/2). The truncated square root of S
    has at least FLINT_BITS n bits (sqrt(S) >= X 2^d); its top n limbs
    with the inexactness of the remainder and the dropped bits give the
    correctly rounded result in every rounding mode. At one and two limbs
    (64-bit), the root of the top limbs of S is first approximated with a
    guard limb, which certifies the truncation in all but rare cases.
*/
/*
    The same with high squarings, for n >= NFLOAT_HYPOT_SQRHIGH_MIN_LIMBS
    (or for exponent gaps d >= 3/8 FLINT_BITS n from
    NFLOAT_HYPOT_SQRHIGH_MIN_LIMBS_GAP limbs):
    with A = the n + 1 top limbs of X^2 (flint_mpn_sqrhigh: a lower bound
    for X^2 / B^(n-1) within 2n - 3 units) and C those of Y^2 shifted
    right by 2d bits (within 2n - 2 units), T = A + C is a lower bound for
    V = (X^2 + Y^2 2^(-2d)) / B^(n-1) within Etot = 4n - 5 units, and
    sqrt(x^2 + y^2) = sqrt(V) B^((n-1)/2) 2^(e_x - FLINT_BITS n). T >= 2^(wn + w - 2)
    roughly (w = FLINT_BITS) is normalized to P = T 2^(2s) of 2n + 2 limbs
    with at least w (2n + 2) - 1 bits, so that R = floor(sqrt(P)) has n + 1
    limbs (top bit set) and differs from sqrt(V 2^(2s)) by at most
    Etot 2^(2s) / (2 sqrt(P)) + 1 <= 2^(3/2) Etot + 1 < 12n units of its
    guard limb. The truncation to n limbs is certified, and inexact,
    unless the guard limb is within 16n + 4 of 0 or B; otherwise
    NFLOAT_HYPOT_DECLINE (the exact computation follows).
*/
/* (measured: with exponents within a few bits of each other, from about
   12 limbs; with a gap of a third of the precision or more, where the
   exact sum S grows to 2.75n limbs or more, from about 5) */
#ifndef NFLOAT_HYPOT_SQRHIGH_MIN_LIMBS
#define NFLOAT_HYPOT_SQRHIGH_MIN_LIMBS 12
#endif
#define NFLOAT_HYPOT_SQRHIGH_MIN_LIMBS_GAP 5
#define NFLOAT_HYPOT_DECLINE (-1)

static int
_nfloat_hypot_sqrhigh(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong d, gr_ctx_t ctx)
{
    ulong H[2 * NFLOAT_MAX_LIMBS];
    ulong A[NFLOAT_MAX_LIMBS + 2], C[NFLOAT_MAX_LIMBS + 2];
    ulong P[2 * NFLOAT_MAX_LIMBS + 2], R[NFLOAT_MAX_LIMBS + 1], M[NFLOAT_MAX_LIMBS];
    slong n = NFLOAT_CTX_NLIMBS(ctx), tn, tl, s2, ls, exp;
    unsigned int bs;
    ulong g, delta;

    /* A, C with the extra limb at index 0 */
    A[0] = flint_mpn_sqrhigh(H, NFLOAT_D(x), n);
    flint_mpn_copyi(A + 1, H, n);
    A[n + 1] = 0;
    C[0] = flint_mpn_sqrhigh(H, NFLOAT_D(y), n);
    flint_mpn_copyi(C + 1, H, n);

    /* C >>= 2d (2d <= FLINT_BITS n + 4 < FLINT_BITS (n + 1)) */
    ls = (2 * d) / FLINT_BITS;
    bs = (2 * d) % FLINT_BITS;
    if (ls != 0)
    {
        flint_mpn_copyi(C, C + ls, n + 1 - ls);
        flint_mpn_zero(C + n + 1 - ls, ls);
    }
    if (bs != 0)
        mpn_rshift(C, C, n + 1, bs);

    /* T = A + C in (A, n + 2) */
    A[n + 1] = mpn_add_n(A, A, C, n + 1);

    tn = (A[n + 1] != 0) ? n + 2 : n + 1;
    tl = FLINT_BITS * tn - flint_clz(A[tn - 1]);

    /* P = T 2^(2s), 2n + 2 limbs, w (2n + 2) - 1 or w (2n + 2) bits */
    s2 = (FLINT_BITS * (2 * n + 2) - tl) & ~WORD(1);
    ls = s2 / FLINT_BITS;
    bs = s2 % FLINT_BITS;
    flint_mpn_zero(P, 2 * n + 2);
    if (bs != 0)
    {
        ulong cy = mpn_lshift(P + ls, A, tn, bs);
        if (ls + tn < 2 * n + 2)
            P[ls + tn] = cy;
    }
    else
        flint_mpn_copyi(P + ls, A, tn);

    flint_mpn_sqrtrem(R, NULL, P, 2 * n + 2);

    g = R[0];
    delta = 16 * n + 4;
    if (g <= delta || g >= UWORD_MAX - delta)
        return NFLOAT_HYPOT_DECLINE;

    flint_mpn_copyi(M, R + 1, n);
    exp = FLINT_BITS + (FLINT_BITS / 2) * (n - 1) - s2 / 2 + NFLOAT_EXP(x);

    if (nfloat_should_round_up(0, ctx))
        NFLOAT_MANT_INCREMENT(M, n, exp);

    _nfloat_copy_limbs(NFLOAT_D(res), M, n);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = 0;
    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

static int
_nfloat_hypot_exact(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong d, gr_ctx_t ctx)
{
    ulong S[3 * NFLOAT_MAX_LIMBS + 4];
    ulong T[2 * NFLOAT_MAX_LIMBS];
    ulong R[2 * NFLOAT_MAX_LIMBS + 4];
    ulong M[NFLOAT_MAX_LIMBS + 1];
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    slong ls = (2 * d) / FLINT_BITS, L, rn, exp, i;
    unsigned int bs = (2 * d) % FLINT_BITS, c;
    int inexact;

#if FLINT_BITS == 64
    if (n == 1)
    {
        /* S = (s3, s2, s1, s0) in registers, then as below */
        ulong h, l, s3, s2, s1, s0, m, rest;
        unsigned int sh = 2 * d, tt;
        slong E;

        umul_ppmm(h, l, NFLOAT_D(x)[0], NFLOAT_D(x)[0]);
        s3 = 0;
        if (sh >= FLINT_BITS)
        {
            /* 2d <= FLINT_BITS + 4 */
            sh -= FLINT_BITS;
            s0 = 0;
            s1 = l << sh;
            s2 = (h << sh) | ((sh == 0) ? 0 : (l >> (FLINT_BITS - sh)));
            s3 = (sh == 0) ? 0 : (h >> (FLINT_BITS - sh));
        }
        else
        {
            s0 = l << sh;
            s1 = (h << sh) | ((sh == 0) ? 0 : (l >> (FLINT_BITS - sh)));
            s2 = (sh == 0) ? 0 : (h >> (FLINT_BITS - sh));
        }
        umul_ppmm(h, l, NFLOAT_D(y)[0], NFLOAT_D(y)[0]);
        add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), UWORD(0), h, l);

        /* the top two limbs, normalized (E = the bit length of S), and
           whether any lower bit of S is nonzero */
        if (s3 != 0)
        {
            c = flint_clz(s3);
            E = 4 * FLINT_BITS - c;
            rest = s0 | ((c == 0) ? s1 : (s1 << c));
            if (c != 0)
            {
                s3 = (s3 << c) | (s2 >> (FLINT_BITS - c));
                s2 = (s2 << c) | (s1 >> (FLINT_BITS - c));
            }
            h = s3; l = s2;
        }
        else if (s2 != 0)
        {
            c = flint_clz(s2);
            E = 3 * FLINT_BITS - c;
            rest = (c == 0) ? s0 : (s0 << c);
            if (c != 0)
            {
                s2 = (s2 << c) | (s1 >> (FLINT_BITS - c));
                s1 = (s1 << c) | (s0 >> (FLINT_BITS - c));
            }
            h = s2; l = s1;
        }
        else
        {
            c = flint_clz(s1);
            E = 2 * FLINT_BITS - c;
            rest = 0;
            if (c != 0)
            {
                s1 = (s1 << c) | (s0 >> (FLINT_BITS - c));
                s0 = s0 << c;
            }
            h = s1; l = s0;
        }

        /* the square root of (h, l) 2^-tt, truncated exactly: the floor of
           the square root of the integer part, which is the floor of the
           square root also with the discarded fraction added */
        tt = E & 1;
        if (tt)
        {
            rest |= l & 1;
            l = (l >> 1) | (h << (FLINT_BITS - 1));
            h = h >> 1;
        }
        inexact = _nfloat_isqrt_2(&m, h, l) | (rest != 0);
        exp = (E + (slong) tt) / 2 + NFLOAT_EXP(x) - FLINT_BITS - d;
        if (inexact && nfloat_should_round_up(0, ctx))
            NFLOAT_MANT_INCREMENT(&m, 1, exp);
        NFLOAT_D(res)[0] = m;
        NFLOAT_EXP(res) = exp;
        NFLOAT_SGNBIT(res) = 0;
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }

    if (n == 2)
    {
        /* S' = X^2 + Y^2 2^-2d = (s4, ..., s0) in registers (the second
           term truncated to the four limbs of the first, an error below
           one unit of s0); U = its top three limbs, normalized, so that
           hypot = sqrt(U 2^-tt) 2^((E + tt) / 2 + EXP(x) - 2 FLINT_BITS)
           for E = the bit length of S'. The root at three limbs is within
           NFLOAT_SMALL_SQRT_ERR_3 units plus at most two for the
           truncations (as below), which certifies the two-limb truncation
           unless the guard limb is within that distance of 0 or B. */
        ulong x3, x2, x1, x0, y3, y2, y1, y0, s4, s3, s2, s1, s0, u2, u1, u0, r2, r1, r0;
        ulong err = NFLOAT_SMALL_SQRT_ERR_3 + 2;
        unsigned int k = 2 * d, b;
        slong E;

        FLINT_MPN_SQR_2X2(x3, x2, x1, x0, NFLOAT_D(x)[1], NFLOAT_D(x)[0]);
        FLINT_MPN_SQR_2X2(y3, y2, y1, y0, NFLOAT_D(y)[1], NFLOAT_D(y)[0]);

        /* (y3, ..., y0) >>= k, k <= 2 FLINT_BITS + 4 */
        if (k >= 2 * FLINT_BITS)
        {
            y0 = y2; y1 = y3; y2 = y3 = 0;
            k -= 2 * FLINT_BITS;
        }
        else if (k >= FLINT_BITS)
        {
            y0 = y1; y1 = y2; y2 = y3; y3 = 0;
            k -= FLINT_BITS;
        }
        b = k;
        if (b != 0)
        {
            y0 = (y0 >> b) | (y1 << (FLINT_BITS - b));
            y1 = (y1 >> b) | (y2 << (FLINT_BITS - b));
            y2 = (y2 >> b) | (y3 << (FLINT_BITS - b));
            y3 = y3 >> b;
        }

        add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, UWORD(0), x3, x2, x1, x0, UWORD(0), y3, y2, y1, y0);

        if (s4 != 0)
        {
            E = 4 * FLINT_BITS + 1;
            u2 = (UWORD(1) << (FLINT_BITS - 1)) | (s3 >> 1);
            u1 = (s3 << (FLINT_BITS - 1)) | (s2 >> 1);
            u0 = (s2 << (FLINT_BITS - 1)) | (s1 >> 1);
        }
        else if (s3 >> (FLINT_BITS - 1))
        {
            E = 4 * FLINT_BITS;
            u2 = s3; u1 = s2; u0 = s1;
        }
        else
        {
            E = 4 * FLINT_BITS - 1;
            u2 = (s3 << 1) | (s2 >> (FLINT_BITS - 1));
            u1 = (s2 << 1) | (s1 >> (FLINT_BITS - 1));
            u0 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
        }

        _nfloat_small_sqrt_3(&r2, &r1, &r0, u2, u1, u0, E & 1);
        if (r0 > err && r0 < UWORD_MAX - err)
        {
            M[0] = r1;
            M[1] = r2;
            exp = (E + (E & 1)) / 2 + NFLOAT_EXP(x) - 2 * FLINT_BITS;
            if (nfloat_should_round_up(0, ctx))
                NFLOAT_MANT_INCREMENT(M, 2, exp);
            NFLOAT_D(res)[0] = M[0];
            NFLOAT_D(res)[1] = M[1];
            NFLOAT_EXP(res) = exp;
            NFLOAT_SGNBIT(res) = 0;
            NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
            return GR_SUCCESS;
        }
    }
#endif

    if (n >= NFLOAT_HYPOT_SQRHIGH_MIN_LIMBS ||
        (n >= NFLOAT_HYPOT_SQRHIGH_MIN_LIMBS_GAP && 8 * d >= 3 * FLINT_BITS * n))
    {
        int status = _nfloat_hypot_sqrhigh(res, x, y, d, ctx);
        if (status != NFLOAT_HYPOT_DECLINE)
            return status;
    }

    L = 2 * n + ls + 1;
    _nfloat_zero_limbs(S, L);
    flint_mpn_sqr(T, NFLOAT_D(x), n);
    if (bs != 0)
        S[ls + 2 * n] = mpn_lshift(S + ls, T, 2 * n, bs);
    else
        _nfloat_copy_limbs(S + ls, T, 2 * n);
    flint_mpn_sqr(T, NFLOAT_D(y), n);
    mpn_add(S, S, L, T, 2 * n);
    while (S[L - 1] == 0)
        L--;

#if FLINT_BITS == 64
    if (n == 2)
    {
        /* U = S 2^-(FLINT_BITS L - c) in [1/2, 1) truncated to n + 1 limbs
           (L >= 2n >= n + 1); sqrt(S) = sqrt(U 2^-t) 2^((FLINT_BITS L - c + t) / 2)
           with t = c mod 2, the root at n + 1 limbs within err units (the
           truncation of U adds at most one): the truncation to n limbs is
           certified, and inexact, unless the guard limb is within err of 0
           or B */
        ulong u2, u1, u0, r2, r1, r0, err;
        ulong t;

        c = flint_clz(S[L - 1]);
        u2 = S[L - 1];
        u1 = S[L - 2];
        u0 = (L >= 3) ? S[L - 3] : 0;
        if (c != 0)
        {
            u2 = (u2 << c) | (u1 >> (FLINT_BITS - c));
            u1 = (u1 << c) | (u0 >> (FLINT_BITS - c));
            u0 = (u0 << c) | ((L >= 4) ? (S[L - 4] >> (FLINT_BITS - c)) : 0);
        }
        t = c & 1;

        _nfloat_small_sqrt_3(&r2, &r1, &r0, u2, u1, u0, t);
        err = NFLOAT_SMALL_SQRT_ERR_3 + 1;
        if (r0 > err && r0 < UWORD_MAX - err)
        {
            M[0] = r1;
            M[1] = r2;
            exp = (FLINT_BITS * L - (slong) c + (slong) t) / 2 + NFLOAT_EXP(x) - FLINT_BITS * n - d;
            if (nfloat_should_round_up(0, ctx))
                NFLOAT_MANT_INCREMENT(M, n, exp);
            _nfloat_copy_limbs(NFLOAT_D(res), M, n);
            NFLOAT_EXP(res) = exp;
            NFLOAT_SGNBIT(res) = 0;
            NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
            return GR_SUCCESS;
        }
    }
#endif

    inexact = (flint_mpn_sqrtrem(R, NULL, S, L) != 0);
    rn = (L + 1) / 2;

    /* the top n limbs of R (normalized) */
    c = flint_clz(R[rn - 1]);
    if (rn == n)
    {
        FLINT_ASSERT(c == 0);
        _nfloat_copy_limbs(M, R, n);
    }
    else
    {
        M[n] = 0;
        if (c != 0)
            mpn_lshift(M, R + rn - n - 1, n + 1, c);
        else
            _nfloat_copy_limbs(M, R + rn - n - 1, n + 1);
        inexact |= (M[0] != 0);
        for (i = 0; i < rn - n - 1 && !inexact; i++)
            inexact = (R[i] != 0);
        flint_mpn_copyi(M, M + 1, n);
    }

    exp = FLINT_BITS * rn - (slong) c + NFLOAT_EXP(x) - FLINT_BITS * n - d;

    if (inexact && nfloat_should_round_up(0, ctx))
        NFLOAT_MANT_INCREMENT(M, n, exp);

    _nfloat_copy_limbs(NFLOAT_D(res), M, n);
    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = 0;
    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
nfloat_hypot(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_abs(res, y, ctx);
        if (NFLOAT_IS_ZERO(y))
            return nfloat_abs(res, x, ctx);
        if (NFLOAT_IS_NAN(x) || NFLOAT_IS_NAN(y))
            return nfloat_nan(res, ctx);
        return nfloat_pos_inf(res, ctx);
    }

    /* |x| (1 + delta), delta <= y^2 / (2 x^2) < 2^(-FLINT_BITS n - 3)
       (in two statements: the operands of | are unsequenced) */
    if (NFLOAT_EXP(x) - NFLOAT_EXP(y) > FLINT_BITS / 2 * NFLOAT_CTX_NLIMBS(ctx) + 2)
    {
        nfloat_abs(res, x, ctx);
        return _nfloat_set_tiny_rel_err(res, res, 1, ctx);
    }
    if (NFLOAT_EXP(y) - NFLOAT_EXP(x) > FLINT_BITS / 2 * NFLOAT_CTX_NLIMBS(ctx) + 2)
    {
        nfloat_abs(res, y, ctx);
        return _nfloat_set_tiny_rel_err(res, res, 1, ctx);
    }

    /* |x| >= |y| */
    if (NFLOAT_EXP(x) > NFLOAT_EXP(y) || (NFLOAT_EXP(x) == NFLOAT_EXP(y) &&
            mpn_cmp(NFLOAT_D(x), NFLOAT_D(y), NFLOAT_CTX_NLIMBS(ctx)) >= 0))
        return _nfloat_hypot_exact(res, x, y, NFLOAT_EXP(x) - NFLOAT_EXP(y), ctx);
    else
        return _nfloat_hypot_exact(res, y, x, NFLOAT_EXP(y) - NFLOAT_EXP(x), ctx);
}
