/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf/elem.h"

/*
    The less common elementary functions of arf (the inverse trigonometric
    and hyperbolic functions and their reciprocal versions, the functions
    divided by pi, sinc of pi x, log10), as compositions of the native
    functions and the arithmetic of arf at a working precision of m limbs,
    retried with more limbs until the rounding is decided (arf/elem.h
    explains the rounding).

    Each composition comes with a bound K for its relative error in units
    of u = 2^(-FLINT_BITS m), counting
      4u for add, sub and mul (the rounding error of ARF_RND_DOWN is below
         an ulp, which is at most 2u; 0 when exact by Sterbenz's lemma, as
         noted),
      2u for div, inv and sqrt (rounded, below an ulp),
      3u for the native functions in ARF_RND_FAST (within an ulp),
      4u for pi (arf_const_pi in ARF_RND_FAST),
    and propagating errors to first order via relative condition numbers
    |x f'(x) / f(x)| (those of the functions used here being at most 1 for
    the arguments where they are applied).  The result R is then rounded
    with an error of K + 2 units of its m-limb mantissa (covering the
    second-order terms).  The special values (zeros, poles, +-1 for the
    inverse functions, infinities) are handled before.

    The composition of an exact dyadic argument is transcendental except
    at the special values (Lindemann-Weierstrass for the inverse
    functions, Niven's theorem for the functions divided by pi), so the
    retries decide the rounding.
*/

typedef enum
{
    ELEM_LOG10,
    ELEM_SINC_PI,
    ELEM_ASIN, ELEM_ACOS, ELEM_ATAN, ELEM_ACOT, ELEM_ASEC, ELEM_ACSC,
    ELEM_ASINH, ELEM_ACOSH, ELEM_ATANH, ELEM_ACOTH, ELEM_ASECH, ELEM_ACSCH
}
elem_func_t;

/* error bounds K (see above), not including division by pi */
static const unsigned char _elem_err[] =
{
    8,
    15,
    11, 11, 3, 3, 11, 11,
    23, 15, 9, 9, 17, 25
};

/* division by pi: 4u for pi, 2u for the division */
#define ELEM_ERR_DIV_PI 6

/* f(x), optionally divided by pi, at m limbs into R (|R| within K units
   of its relative error, see above) */
static void
_arf_elem_inner(arf_t R, elem_func_t f, int div_pi, const arf_t x, slong m)
{
    arf_t X, A, B, S, O;
    slong p = FLINT_BITS * m;
    int sgn = ARF_SGNBIT(x);

    arf_init(X);
    arf_init(A);
    arf_init(B);
    arf_init(S);
    arf_init(O);

    /* X = |x| */
    arf_abs(X, x);
    arf_one(O);

    switch (f)
    {
        case ELEM_LOG10:
            arf_log(A, X, p, ARF_RND_FAST);
            arf_set_ui(B, 10);
            arf_log(B, B, p, ARF_RND_FAST);
            arf_div(R, A, B, p, ARF_RND_DOWN);
            break;

        case ELEM_SINC_PI:
            /* sin(pi x) / (pi x): 3u + 4u + 4u + 2u */
            arf_sin_pi(A, X, p, ARF_RND_FAST);
            arf_const_pi(B, p, ARF_RND_FAST);
            arf_mul(B, B, X, p, ARF_RND_DOWN);
            arf_div(R, A, B, p, ARF_RND_DOWN);
            break;

        case ELEM_ASIN:
        case ELEM_ACOS:
        case ELEM_ASEC:
        case ELEM_ACSC:
            /* S = sqrt((1 - |x|)(1 + |x|)) or sqrt((|x| - 1)(|x| + 1)),
               8u (the difference is exact for 1/2 <= |x| <= 2) */
            if (f == ELEM_ASIN || f == ELEM_ACOS)
                arf_sub(A, O, X, p, ARF_RND_DOWN);
            else
                arf_sub(A, X, O, p, ARF_RND_DOWN);
            arf_add(B, O, X, p, ARF_RND_DOWN);
            arf_mul(A, A, B, p, ARF_RND_DOWN);
            arf_sqrt(S, A, p, ARF_RND_DOWN);

            /* asin x = atan2(x, S), acos x = atan2(S, x),
               asec x = atan2(S, sgn(x)), acsc x = atan2(sgn(x), S) */
            if (f == ELEM_ASEC || f == ELEM_ACSC)
                arf_set_si(X, sgn ? -1 : 1);
            else
                arf_set(X, x);
            if (f == ELEM_ASIN || f == ELEM_ACSC)
                arf_atan2(R, X, S, p, ARF_RND_FAST);
            else
                arf_atan2(R, S, X, p, ARF_RND_FAST);
            break;

        case ELEM_ATAN:
            arf_atan(R, x, p, ARF_RND_FAST);
            break;

        case ELEM_ACOT:
            /* acot x = atan2(sgn(x), |x|) */
            arf_set_si(A, sgn ? -1 : 1);
            arf_atan2(R, A, X, p, ARF_RND_FAST);
            break;

        case ELEM_ACSCH:
            /* asinh(1/|x|), 2u + 23u */
            arf_ui_div(X, 1, X, p, ARF_RND_DOWN);
            FLINT_FALLTHROUGH;

        case ELEM_ASINH:
            /* log1p(y + y^2 / (1 + sqrt(1 + y^2))), y = |x| */
            arf_mul(A, X, X, p, ARF_RND_DOWN);                /* 4u */
            arf_add(B, A, O, p, ARF_RND_DOWN);                /* 8u */
            arf_sqrt(B, B, p, ARF_RND_DOWN);                  /* 6u */
            arf_add(B, B, O, p, ARF_RND_DOWN);                /* 10u */
            arf_div(A, A, B, p, ARF_RND_DOWN);                /* 16u */
            arf_add(A, A, X, p, ARF_RND_DOWN);                /* 20u */
            arf_log1p(R, A, p, ARF_RND_FAST);                 /* 23u */
            if (sgn)
                arf_neg(R, R);
            break;

        case ELEM_ACOSH:
            /* log1p(d + sqrt(d (x + 1))), d = x - 1 */
            arf_sub(A, X, O, p, ARF_RND_DOWN);                /* 4u, exact for x <= 2 */
            arf_add(B, X, O, p, ARF_RND_DOWN);                /* 4u */
            arf_mul(B, A, B, p, ARF_RND_DOWN);                /* 12u */
            arf_sqrt(B, B, p, ARF_RND_DOWN);                  /* 8u */
            arf_add(A, A, B, p, ARF_RND_DOWN);                /* 12u */
            arf_log1p(R, A, p, ARF_RND_FAST);                 /* 15u */
            break;

        case ELEM_ATANH:
        case ELEM_ACOTH:
            /* log1p(2|x| / (1 - |x|)) / 2, log1p(2 / (|x| - 1)) / 2 */
            if (f == ELEM_ATANH)
            {
                arf_sub(A, O, X, p, ARF_RND_DOWN);            /* 4u, exact for |x| >= 1/2 */
                arf_div(A, X, A, p, ARF_RND_DOWN);            /* 6u */
            }
            else
            {
                arf_sub(A, X, O, p, ARF_RND_DOWN);            /* 4u, exact for |x| <= 2 */
                arf_ui_div(A, 1, A, p, ARF_RND_DOWN);         /* 6u */
            }
            arf_mul_2exp_si(A, A, 1);
            arf_log1p(R, A, p, ARF_RND_FAST);                 /* 9u */
            arf_mul_2exp_si(R, R, -1);
            if (sgn)
                arf_neg(R, R);
            break;

        case ELEM_ASECH:
            /* log1p((a + sqrt(a (1 + x))) / x), a = 1 - x */
            arf_sub(A, O, X, p, ARF_RND_DOWN);                /* 4u, exact for x >= 1/2 */
            arf_add(B, X, O, p, ARF_RND_DOWN);                /* 4u */
            arf_mul(B, A, B, p, ARF_RND_DOWN);                /* 12u */
            arf_sqrt(B, B, p, ARF_RND_DOWN);                  /* 8u */
            arf_add(A, A, B, p, ARF_RND_DOWN);                /* 12u */
            arf_div(A, A, X, p, ARF_RND_DOWN);                /* 14u */
            arf_log1p(R, A, p, ARF_RND_FAST);                 /* 17u */
            break;
    }

    if (div_pi)
    {
        arf_const_pi(A, p, ARF_RND_FAST);
        arf_div(R, R, A, p, ARF_RND_DOWN);
    }

    arf_clear(X);
    arf_clear(A);
    arf_clear(B);
    arf_clear(S);
    arf_clear(O);
}

/* f(x) (optionally divided by pi) for a finite nonzero x that passed
   the special-value checks of the functions below; res may alias x */
static int
_arf_elem(arf_t res, elem_func_t f, int div_pi, const arf_t x, slong prec, arf_rnd_t rnd)
{
    arf_t R;
    slong extra, m;
    ulong K;
    int status;

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf elementary function: exact precision not supported\n");

#if FLINT_BITS == 64
    /* the inverse functions at one and two limbs by the in-register
       compositions (elem_arc_fast.c), a shorter mantissa padded */
    if (f >= ELEM_ASIN && f <= ELEM_ACSCH && prec <= 2 * FLINT_BITS
        && !COEFF_IS_MPZ(ARF_EXP(x)) && ARF_SIZE(x) <= 2)
    {
        static const unsigned char arc[] = {
            ARF_ARC_ASIN, ARF_ARC_ACOS, ARF_ARC_ATAN, ARF_ARC_ACOT,
            ARF_ARC_ASEC, ARF_ARC_ACSC, ARF_ARC_ASINH, ARF_ARC_ACOSH,
            ARF_ARC_ATANH, ARF_ARC_ACOTH, ARF_ARC_ASECH, ARF_ARC_ACSCH };
        _arf_arc_out_struct out;
        nn_srcptr d;
        slong xn, n, w;
        ulong dd[2], y[5];

        ARF_GET_MPN_READONLY(d, xn, x);
        n = (xn == 1 && prec <= FLINT_BITS) ? 1 : 2;
        w = n + 1;
        dd[0] = (xn == n) ? d[0] : 0;
        dd[1] = d[xn - 1];

        if (_arf_arc_fast(&out, dd, n, ARF_EXP(x), ARF_SGNBIT(x), arc[f - ELEM_ASIN], div_pi) != ARF_ELEM_UNCERTAIN)
        {
            flint_mpn_copyi(y, out.T, w + 1);
            status = _arf_set_mpn_err(res, y, w + 1, out.k - FLINT_BITS * w, ARF_ELEM_EB(out.err), ARF_ELEM_EB(out.err), out.sgnbit, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
#endif

    K = _elem_err[f] + (div_pi ? ELEM_ERR_DIV_PI : 0);
    arf_init(R);

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_srcptr rd;
        nn_ptr y;
        slong rn;
        TMP_INIT;

        /* m limbs: the working precision plus the bits of the error */
        m = _arf_elem_limbs(prec, rnd, extra + FLINT_BIT_COUNT(K + 2));

        _arf_elem_inner(R, f, div_pi, x, m);

        if (arf_is_special(R))
        {
            /* (not reached for the arguments the callers pass) */
            arf_nan(res);
            status = 0;
            break;
        }

        /* the mantissa padded to m limbs: within K + 2 units */
        ARF_GET_MPN_READONLY(rd, rn, R);
        TMP_START;
        y = TMP_ALLOC((m + 1) * sizeof(ulong));
        _flint_mpn_zero_small(y, m - rn);
        flint_mpn_copyi(y + m - rn, rd, rn);
        if (!COEFF_IS_MPZ(ARF_EXP(R)))
            status = _arf_set_mpn_err(res, y, m, ARF_EXP(R) - FLINT_BITS * m, ARF_ELEM_EB(K + 2), ARF_ELEM_EB(K + 2), ARF_SGNBIT(R), prec, rnd);
        else
        {
            /* a bignum exponent: rounded with a zero unit exponent, then
               shifted (res may alias x, which R does not) */
            fmpz_t e;
            fmpz_init(e);
            fmpz_sub_ui(e, ARF_EXPREF(R), FLINT_BITS * m);
            status = _arf_set_mpn_err(res, y, m, 0, ARF_ELEM_EB(K + 2), ARF_ELEM_EB(K + 2), ARF_SGNBIT(R), prec, rnd);
            if (status == 1)
                fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), e);
            fmpz_clear(e);
        }
        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            break;
    }

    arf_clear(R);
    return status;
}

/* The public functions: the special values, the arguments where the
   value is a dyadic number times (1 + delta) with a tiny delta of known
   sign (which the compositions could never round, their error bounds
   being two-sided), then the compositions. */

/* the exponent of x, a bignum clamped beyond every threshold below */
static slong
_exp_clamped(const arf_t x, slong prec)
{
    slong lim = prec + FLINT_BITS * (ARF_SIZE(x) + 2) + 64;
    if (!COEFF_IS_MPZ(ARF_EXP(x)))
        return ARF_EXP(x);
    return (fmpz_sgn(ARF_EXPREF(x)) < 0) ? -lim : lim;
}

/* For x = 1 -+ eps (sgn = -1 resp. 1) with 0 < eps <= 2^-prec: the
   functions acos, acosh, asec and asech are sqrt(2 eps)(1 + delta) with
   0 < |delta| < 2^-P, P = 3 - e (acos: delta > 0, acosh: delta < 0) resp.
   P = 1 - e (asech: delta > 0, asec: delta < 0) for eps < 2^e.  When
   sqrt(2 eps) is a dyadic number (the exact root at the length of eps),
   the result is rounded from it; otherwise, or when the rounding is not
   decided, ARF_ELEM_UNCERTAIN (the compositions decide, sqrt(2 eps)
   being irrational then). */
static int
_arf_near_one(arf_t res, const arf_t x, int sgn, int side, int p3, slong prec, arf_rnd_t rnd)
{
    arf_t eps, s;
    slong e;
    int status = ARF_ELEM_UNCERTAIN;

    /* (eps < 1/2 needs x in [1/2, 1) resp. (1, 2): the exponent 0 resp. 1) */
    if (COEFF_IS_MPZ(ARF_EXP(x)) || ARF_EXP(x) != ((sgn < 0) ? 0 : 1))
        return ARF_ELEM_UNCERTAIN;

    /* eps <= 2^-prec needs the top min(prec, FLINT_BITS) bits of the
       mantissa all ones (x < 1) resp. a leading one followed by zeros
       (x > 1): a cheap test before the exact difference */
    {
        ulong top = ARF_NOPTR_D(x)[0];
        slong b = FLINT_MIN(prec, FLINT_BITS - 1);
        if (ARF_SIZE(x) > ARF_NOPTR_LIMBS)
            top = ARF_PTR_D(x)[ARF_SIZE(x) - 1];
        else
            top = ARF_NOPTR_D(x)[ARF_SIZE(x) - 1];
        if (sgn < 0 ? ((top >> (FLINT_BITS - b)) != (UWORD_MAX >> (FLINT_BITS - b)))
                    : (((top << 1) >> (FLINT_BITS - b)) != 0))
            return ARF_ELEM_UNCERTAIN;
    }

    arf_init(eps);
    arf_init(s);

    /* eps = |x - 1|, exact */
    arf_sub_ui(eps, x, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_abs(eps, eps);
    FLINT_ASSERT(!arf_is_zero(eps));
    e = _exp_clamped(eps, prec);

    /* (P >= prec + 3, as _arf_set_tiny_rel_err needs: e <= -prec for
       P = 3 - e, e <= -prec - 2 for P = 1 - e) */
    if (e <= -prec - (p3 ? 0 : 2))
    {
        arf_mul_2exp_si(eps, eps, 1);
        if (arf_sqrt(s, eps, FLINT_BITS * (ARF_SIZE(eps) + 1), ARF_RND_DOWN) == 0)
            status = _arf_set_tiny_rel_err(res, s, side, (p3 ? 3 : 1) - e, prec, rnd);
    }

    arf_clear(eps);
    arf_clear(s);
    return status;
}


int
arf_log10(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x) || arf_sgn(x) < 0)
    {
        if (arf_is_zero(x))
            arf_neg_inf(res);
        else if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return 0;
    }

    /* x = 10^k exactly: m 2^s with m = 5^s (2 s < bits(m) <= 3 s) */
    {
        fmpz_t m, s;
        int exact = 0;
        slong k = 0;

        fmpz_init(m);
        fmpz_init(s);
        arf_get_fmpz_2exp(m, s, x);
        if (fmpz_sgn(s) >= 0 && fmpz_fits_si(s) && fmpz_get_si(s) <= WORD(1) << 20)
        {
            slong b = fmpz_bits(m);
            k = fmpz_get_si(s);
            if (k == 0)
                exact = fmpz_is_one(m);
            else if (2 * k < b && b <= 3 * k)
            {
                fmpz_ui_pow_ui(s, 5, k);
                exact = fmpz_equal(m, s);
            }
        }
        fmpz_clear(m);
        fmpz_clear(s);

        if (exact)
            return arf_set_round_si(res, k, prec, arf_rnd_relaxed_to_strict(rnd));
    }

    return _arf_elem(res, ELEM_LOG10, 0, x, prec, rnd);
}

int
arf_sinc_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    slong e;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_one(res);
        else
            arf_nan(res);
        return 0;
    }

    /* zero at the nonzero integers */
    if (arf_is_int(x))
    {
        arf_zero(res);
        return 0;
    }

    /* 1 - delta, delta < (pi x)^2 / 6 < 2 x^2 */
    e = _exp_clamped(x, prec);
    if (2 * e <= -prec - 2)
        return _arf_one_plus_tiny(res, 0, 1, prec, rnd);

    return _arf_elem(res, ELEM_SINC_PI, 0, x, prec, rnd);
}

/* the inverse trigonometric functions, with the division by pi */

static int
_arf_asin(arf_t res, const arf_t x, int div_pi, slong prec, arf_rnd_t rnd)
{
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else
            arf_nan(res);
        return 0;
    }

    c = arf_cmpabs_ui(x, 1);
    if (c > 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
    {
        /* +-pi/2 */
        if (div_pi)
        {
            arf_set_si(res, arf_sgn(x));
            arf_mul_2exp_si(res, res, -1);
            return 0;
        }
        return arf_const_pi_mul_2exp(res, -1, arf_sgn(x) < 0, prec, rnd);
    }

    if (!div_pi)
    {
        /* asin x = x (1 + delta), 0 < delta < x^2 / 5 < 2^(2e - 2) */
        slong e = _exp_clamped(x, prec);
        if (2 * e <= -prec - 2)
        {
            int status = _arf_set_tiny_rel_err(res, x, 1, 2 - 2 * e, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }
    else if (!COEFF_IS_MPZ(ARF_EXP(x)) && ARF_EXP(x) == 0)
    {
        /* asin(+-(1 - eps)) / pi = +-(1/2)(1 - delta) with 0 < delta =
           2 acos(1 - eps) / pi < sqrt(2 eps) < 2^((e + 1) / 2 + 1) for
           eps < 2^e */
        arf_t eps;
        slong e;
        arf_init(eps);
        arf_sub_ui(eps, x, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_abs(eps, eps);
        e = _exp_clamped(eps, prec);
        arf_clear(eps);
        if (-(e + 1) / 2 - 1 >= prec + 3)
        {
            arf_t h;
            int status;
            arf_init(h);
            arf_set_ui_2exp_si(h, 1, -1);
            if (arf_sgn(x) < 0)
                arf_neg(h, h);
            status = _arf_set_tiny_rel_err(res, h, -1, -(e + 1) / 2 - 1, prec, rnd);
            arf_clear(h);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    return _arf_elem(res, ELEM_ASIN, div_pi, x, prec, rnd);
}

static int
_arf_acos(arf_t res, const arf_t x, int div_pi, slong prec, arf_rnd_t rnd)
{
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            /* pi/2 */
            if (div_pi)
            {
                arf_set_ui_2exp_si(res, 1, -1);
                return 0;
            }
            return arf_const_pi_mul_2exp(res, -1, 0, prec, rnd);
        }
        arf_nan(res);
        return 0;
    }

    c = arf_cmpabs_ui(x, 1);
    if (c > 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
    {
        /* 0 or pi */
        if (arf_sgn(x) > 0)
        {
            arf_zero(res);
            return 0;
        }
        if (div_pi)
        {
            arf_one(res);
            return 0;
        }
        return arf_const_pi_mul_2exp(res, 0, 0, prec, rnd);
    }

    if (!div_pi && arf_sgn(x) > 0)
    {
        /* acos(1 - eps) = sqrt(2 eps)(1 + delta), 0 < delta < eps / 8 */
        int status = _arf_near_one(res, x, -1, 1, 1, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    if (div_pi && arf_sgn(x) < 0 && !COEFF_IS_MPZ(ARF_EXP(x)) && ARF_EXP(x) == 0)
    {
        /* acos(-1 + eps) / pi = 1 - delta with 0 < delta = acos(1 - eps) /
           pi < sqrt(2 eps) / 2 < 2^((e + 1) / 2) for eps < 2^e */
        arf_t eps;
        slong e;
        arf_init(eps);
        arf_add_ui(eps, x, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
        e = _exp_clamped(eps, prec);
        arf_clear(eps);
        if ((e + 1) / 2 <= -prec - 1)
            return _arf_one_plus_tiny(res, 0, 1, prec, rnd);
    }

    if (div_pi)
    {
        /* acos(x) / pi = (1/2)(1 - 2 asin(x) / pi) with |2 asin(x) / pi| <
           2^(e + 1) for a tiny x */
        slong e = _exp_clamped(x, prec);
        if (e <= -prec - 4)
        {
            arf_t h;
            int status;
            arf_init(h);
            arf_set_ui_2exp_si(h, 1, -1);
            status = _arf_set_tiny_rel_err(res, h, (arf_sgn(x) > 0) ? -1 : 1, -e - 1, prec, rnd);
            arf_clear(h);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    return _arf_elem(res, ELEM_ACOS, div_pi, x, prec, rnd);
}

int
arf_asin(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_asin(res, x, 0, prec, rnd);
}

int
arf_acos(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_acos(res, x, 0, prec, rnd);
}

int
arf_asin_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_asin(res, x, 1, prec, rnd);
}

int
arf_acos_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    return _arf_acos(res, x, 1, prec, rnd);
}

int
arf_atan_pi(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else if (arf_is_nan(x))
            arf_nan(res);
        else
        {
            /* +-1/2 */
            arf_set_si(res, arf_sgn(x));
            arf_mul_2exp_si(res, res, -1);
        }
        return 0;
    }

    /* +-1/4 at +-1 */
    if (arf_cmpabs_ui(x, 1) == 0)
    {
        arf_set_si(res, arf_sgn(x));
        arf_mul_2exp_si(res, res, -2);
        return 0;
    }

    /* atan(x) / pi = (1/2)(1 - delta) for x > 0 with
       0 < delta = 2 atan(1/x) / pi < 2^(1 - e) */
    {
        slong e = _exp_clamped(x, prec);
        if (e - 1 >= prec + 3)
        {
            arf_t h;
            int status;
            arf_init(h);
            arf_set_ui_2exp_si(h, 1, -1);
            if (arf_sgn(x) < 0)
                arf_neg(h, h);
            status = _arf_set_tiny_rel_err(res, h, -1, e - 1, prec, rnd);
            arf_clear(h);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    return _arf_elem(res, ELEM_ATAN, 1, x, prec, rnd);
}

int
arf_acot(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return arf_const_pi_mul_2exp(res, -1, 0, prec, rnd);
        if (arf_is_nan(x))
            arf_nan(res);
        else
            arf_zero(res);
        return 0;
    }

    /* acot x = (1/x)(1 - delta), 0 < delta < 1 / (3 x^2) < 2^(1 - 2e) */
    {
        slong e = _exp_clamped(x, prec);
        if (2 * e - 1 >= prec + 3)
        {
            int status = _arf_set_tiny_recip_rel_err(res, x, -1, 2 * e - 1, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    return _arf_elem(res, ELEM_ACOT, 0, x, prec, rnd);
}

int
arf_asec(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_inf(x))
            return arf_const_pi_mul_2exp(res, -1, 0, prec, rnd);
        arf_nan(res);
        return 0;
    }

    c = arf_cmpabs_ui(x, 1);
    if (c < 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
    {
        if (arf_sgn(x) > 0)
        {
            arf_zero(res);
            return 0;
        }
        return arf_const_pi_mul_2exp(res, 0, 0, prec, rnd);
    }

    if (arf_sgn(x) > 0)
    {
        /* asec(1 + eps) = sqrt(2 eps)(1 - delta), 0 < delta < eps / 2 */
        int status = _arf_near_one(res, x, 1, -1, 0, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    return _arf_elem(res, ELEM_ASEC, 0, x, prec, rnd);
}

int
arf_acsc(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_inf(x))
            arf_zero(res);
        else
            arf_nan(res);
        return 0;
    }

    c = arf_cmpabs_ui(x, 1);
    if (c < 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
        return arf_const_pi_mul_2exp(res, -1, arf_sgn(x) < 0, prec, rnd);

    /* acsc x = (1/x)(1 + delta), 0 < delta < 1 / (5 x^2) < 2^(1 - 2e) */
    {
        slong e = _exp_clamped(x, prec);
        if (2 * e - 1 >= prec + 3)
        {
            int status = _arf_set_tiny_recip_rel_err(res, x, 1, 2 * e - 1, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    return _arf_elem(res, ELEM_ACSC, 0, x, prec, rnd);
}

/* the inverse hyperbolic functions */

int
arf_asinh(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    slong e;

    if (arf_is_special(x))
    {
        arf_set(res, x);
        return 0;
    }

    /* asinh x = x (1 - delta), 0 < delta < x^2 / 6 < 2^(2e - 2) */
    e = _exp_clamped(x, prec);
    if (2 * e <= -prec - 2)
    {
        int status = _arf_set_tiny_rel_err(res, x, -1, 2 - 2 * e, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    return _arf_elem(res, ELEM_ASINH, 0, x, prec, rnd);
}

int
arf_acosh(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return 0;
    }

    c = arf_cmp_ui(x, 1);
    if (c < 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
    {
        arf_zero(res);
        return 0;
    }

    /* acosh(1 + eps) = sqrt(2 eps)(1 - delta), 0 < delta < eps / 8 */
    {
        int status = _arf_near_one(res, x, 1, -1, 1, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    return _arf_elem(res, ELEM_ACOSH, 0, x, prec, rnd);
}

int
arf_atanh(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    slong e;
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else
            arf_nan(res);
        return 0;
    }

    c = arf_cmpabs_ui(x, 1);
    if (c > 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
    {
        if (arf_sgn(x) > 0)
            arf_pos_inf(res);
        else
            arf_neg_inf(res);
        return 0;
    }

    /* atanh x = x (1 + delta), 0 < delta < x^2 / 2 < 2^(2e - 1) */
    e = _exp_clamped(x, prec);
    if (2 * e <= -prec - 2)
    {
        int status = _arf_set_tiny_rel_err(res, x, 1, 1 - 2 * e, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    return _arf_elem(res, ELEM_ATANH, 0, x, prec, rnd);
}

int
arf_acoth(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_inf(x))
            arf_zero(res);
        else
            arf_nan(res);
        return 0;
    }

    c = arf_cmpabs_ui(x, 1);
    if (c < 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
    {
        if (arf_sgn(x) > 0)
            arf_pos_inf(res);
        else
            arf_neg_inf(res);
        return 0;
    }

    /* acoth x = (1/x)(1 + delta), 0 < delta < 1 / (2 x^2) < 2^(1 - 2e) */
    {
        slong e = _exp_clamped(x, prec);
        if (2 * e - 1 >= prec + 3)
        {
            int status = _arf_set_tiny_recip_rel_err(res, x, 1, 2 * e - 1, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    return _arf_elem(res, ELEM_ACOTH, 0, x, prec, rnd);
}

int
arf_asech(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    int c;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return 0;
    }

    c = arf_cmp_ui(x, 1);
    if (c > 0 || arf_sgn(x) < 0)
    {
        arf_nan(res);
        return 0;
    }
    if (c == 0)
    {
        arf_zero(res);
        return 0;
    }

    /* asech(1 - eps) = sqrt(2 eps)(1 + delta), 0 < delta < eps / 2 */
    {
        int status = _arf_near_one(res, x, -1, 1, 0, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    return _arf_elem(res, ELEM_ASECH, 0, x, prec, rnd);
}

int
arf_acsch(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        if (arf_is_inf(x))
            arf_zero(res);
        else
            arf_nan(res);
        return 0;
    }

    /* acsch x = (1/x)(1 - delta), 0 < delta < 1 / (6 x^2) < 2^(1 - 2e) */
    {
        slong e = _exp_clamped(x, prec);
        if (2 * e - 1 >= prec + 3)
        {
            int status = _arf_set_tiny_recip_rel_err(res, x, -1, 2 * e - 1, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN)
                return status;
        }
    }

    return _arf_elem(res, ELEM_ACSCH, 0, x, prec, rnd);
}
