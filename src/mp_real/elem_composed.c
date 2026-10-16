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
#include "mp_real.h"
#include "impl.h"

/* log1p, expm1, the hyperbolic functions and the inverse trigonometric
   and hyperbolic functions of a ball to a relative accuracy of about
   2^-prec, as compositions of exp, log and atan in ball arithmetic.

   The midpoint m is evaluated exactly as given, by formulas without
   cancellation:

     log1p x = log(1 + x),  1 + x exact or to prec + z bits for |x| < 2^-z,
     expm1 x = exp(x) - 1,  exp(x) to prec + z bits for |x| < 2^-z,
     sinh x = (expm1(x) - expm1(-x)) / 2,  cosh x = 1 + (expm1(x) + expm1(-x)) / 2,
     tanh |x| = expm1(2|x|) / (expm1(2|x|) + 2),
     asin x = atan(x / sqrt((1 - x)(1 + x))),  acos x = 2 atan(sqrt((1 - x) / (1 + x))),
     asinh |x| = log1p(|x| + x^2 / (1 + sqrt(1 + x^2))),
     acosh x = log1p(d + sqrt(d (x + 1))),  d = x - 1,
     atanh |x| = log1p(2|x| / (1 - |x|)) / 2,

   with the differences 1 - |x| and x - 1 exact, and the small and large
   arguments by the leading terms and bounds for the rest (x, log(2|x|)
   resp. +-1).  The radius rho of x then enters as the image of the
   derivative bound over the ball (rho itself for tanh and asinh, rho /
   sqrt(1 - |m| - rho), rho / sqrt(m - 1 - rho), rho / (1 - |m| - rho)
   for asin, acos, acosh, atanh), or through exp and log for the others.
   Outputs may alias the input. */

/* 2^(E - 1) <= |m| < 2^E for a nonzero midpoint */
FLINT_FORCE_INLINE slong
_ec_emid(const mp_real_t m)
{
    return FLINT_BITS * (m->exp - 1) + FLINT_BIT_COUNT(m->d[m->size - 1]);
}

/* the working limbs */
FLINT_FORCE_INLINE slong
_ec_wp(slong prec)
{
    return mp_real_prec_bits(prec + 16);
}

/* limbs for which 1 +- m, m - 1 are exact (|m| >= 1/B), else wp */
FLINT_FORCE_INLINE slong
_ec_exact_limbs(const mp_real_t m, slong wp)
{
    if (m->size != 0 && m->exp >= 0)
        return FLINT_MAX(wp, m->size + 2);
    return wp;
}

/* 1 = |m| (exact): 0 resp. +-1 for |m| < 1, > 1 */
static int
_ec_cmpabs_one(const mp_real_t m)
{
    if (m->size == 0 || m->exp <= 0)
        return -1;
    if (m->exp >= 2)
        return 1;
    /* exp = 1: |m| in [1, B) */
    if (m->d[m->size - 1] > 1)
        return 1;
    return (m->size == 1) ? 0 : 1;
}

/* the radius bound rho / sqrt(D) (root = 1) resp. rho / D (root = 0) for
   D = a - rho, a exact: added to res's radius; returns 0 if D cannot be
   certified positive.  With alt = 1, the alternative bound 2 sqrt(a +
   rho) is used when smaller (for a ball reaching close to the branch
   point, where the function varies like the square root), with alt = 2
   only if a + rho < 1. */
static int
_ec_add_rad_div(mp_real_t res, const mp_real_t a, ulong xerr, slong xanc,
    int root, int alt, slong wp)
{
    mp_real_t R, D, V;
    int ok = 1;

    mp_real_init(R);
    mp_real_init(D);
    mp_real_init(V);

    _mp_real_set_mpn_2exp(R, &xerr, 1, FLINT_BITS * xanc);
    mp_real_sub(D, a, R, wp + 2);

    if (D->size == 0 || D->negative || mp_real_rel_radius_lt_2exp_si(D) > -32)
        ok = 0;
    else
    {
        if (root)
            mp_real_sqrt(D, D, 2);
        mp_real_div(V, R, D, 2);

        if (alt)
        {
            mp_real_add(D, a, R, 2);
            if (alt == 1 || mp_real_abs_bound_lt_2exp_si(D) <= 0)
            {
                mp_real_sqrt(D, D, 2);
                mp_real_mul_2exp_si(D, D, 1);
                if (mp_real_abs_bound_lt_2exp_si(D) < mp_real_abs_bound_lt_2exp_si(V))
                    mp_real_swap(V, D);
            }
        }

        _mp_real_elem_add_mag(res, V);
    }

    mp_real_clear(R);
    mp_real_clear(D);
    mp_real_clear(V);
    return ok;
}

/* the midpoint of x (exact, without the zero padding limbs an inexact
   ball may carry; its limbs are x's) and the radius */
#define EC_SPLIT(mid, x, xerr, xanc) \
    do { \
        mid = *(x); \
        mid.err = 0; \
        while (mid.size > 0 && mid.d[0] == 0) \
        { \
            mid.d++; \
            mid.size--; \
        } \
        if (mid.size == 0) \
            mid.exp = 0; \
        xerr = (x)->err; \
        xanc = (x)->exp - (x)->size; \
    } while (0)

/* log1p *********************************************************************/

int
mp_real_log1p_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_t u, one;
    slong e, z, n;
    int ok;

    prec = FLINT_MAX(prec, 2);

    if (x->size == 0 && x->err == 0)
    {
        mp_real_zero(res);
        return 1;
    }

    /* |y| < 2^e over the ball */
    e = mp_real_abs_bound_lt_2exp_si(x);

    if (e <= -(prec + 4))
    {
        /* |log1p(y) - y| <= y^2 for |y| <= 1/2 */
        mp_real_set(res, x);
        mp_real_add_error_2exp_si(res, _mp_real_err_exp_clamp(2 * e, e + FLINT_BITS, prec));
        return 1;
    }

    /* u = 1 + x, exact or (for |x| < 2^-z) to prec + z + 16 bits: an
       absolute error far below |x| 2^-prec */
    z = (e < 0) ? -e : 0;
    n = mp_real_prec_bits(prec + z + 16);
    if (x->size != 0)
        n = _ec_exact_limbs(x, n);

    mp_real_init(u);
    mp_real_init(one);
    mp_real_set_ui(one, 1);
    mp_real_add(u, x, one, n);
    ok = mp_real_log_bits(res, u, prec);
    mp_real_clear(u);
    mp_real_clear(one);
    return ok;
}

/* expm1 *********************************************************************/

void
mp_real_expm1_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_t E, one;
    slong e;

    prec = FLINT_MAX(prec, 2);

    if (x->size == 0 && x->err == 0)
    {
        mp_real_zero(res);
        return;
    }

    e = mp_real_abs_bound_lt_2exp_si(x);

    if (e <= -(prec + 4))
    {
        /* |expm1(y) - y| <= y^2 for |y| <= 1/2 */
        mp_real_set(res, x);
        mp_real_add_error_2exp_si(res, _mp_real_err_exp_clamp(2 * e, e + FLINT_BITS, prec));
        return;
    }

    /* x <= -2^(FLINT_BITS - 7) over the ball: -1 + exp(x), 0 < exp(x)
       < 2^(-2^(FLINT_BITS - 7)) (exp itself would throw) */
    if (x->size != 0 && x->negative && _ec_emid(x) >= FLINT_BITS - 5
        && mp_real_rel_radius_lt_2exp_si(x) <= -2)
    {
        _mp_real_elem_set_error(res, 1, _mp_real_err_exp_clamp(-(WORD(1) << (FLINT_BITS - 7)), FLINT_BITS, prec));
        mp_real_neg(res, res);
        return;
    }

    mp_real_init(E);
    mp_real_init(one);
    mp_real_set_ui(one, 1);

    if (e <= 0)
    {
        /* |x| < 2^-z, z = -e: exp(x) to prec + z bits */
        slong z = -e;
        mp_real_exp_bits(E, x, prec + z + 8);
        mp_real_sub(res, E, one, mp_real_prec_bits(prec + z + 16));
    }
    else
    {
        /* exp(x) - 1 >= e - 1 resp. in [-1, 1/e - 1]: no cancellation
           (except for a wide ball, which the radius shows) */
        mp_real_exp_bits(E, x, prec + 4);

        /* exp(x) negligible against 1 (large negative x): -1 with the
           bound for exp(x) as the radius, clamped (adding E's tiny
           radius would pad the mantissa down to its scale) */
        if (mp_real_abs_bound_lt_2exp_si(E) < -prec - 2 * FLINT_BITS)
        {
            _mp_real_elem_set_error(res, 1, _mp_real_err_exp_clamp(
                mp_real_abs_bound_lt_2exp_si(E), 1, prec));
            mp_real_neg(res, res);
        }
        else
            mp_real_sub(res, E, one, _ec_wp(prec));
    }

    mp_real_clear(E);
    mp_real_clear(one);
}

/* sinh and cosh ****************************************************************/

void
mp_real_sinh_cosh_bits(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    mp_real_t A, Bn, S, C, one;
    ulong xerr;
    slong xanc, wp;

    prec = FLINT_MAX(prec, 2);
    wp = _ec_wp(prec);
    EC_SPLIT(mid, x, xerr, xanc);

    mp_real_init(S);
    mp_real_init(C);

    /* at the midpoint */
    if (mid.size == 0)
    {
        mp_real_zero(S);
        mp_real_set_ui(C, 1);
    }
    else
    {
        slong e = _ec_emid(&mid);

        if (e <= -(prec / 2 + 4))
        {
            /* |sinh m - m| <= |m|^3, |cosh m - 1| <= m^2 for |m| <= 1/2 */
            mp_real_set(S, &mid);
            mp_real_add_error_2exp_si(S, _mp_real_err_exp_clamp(3 * e, e + FLINT_BITS, prec));
            _mp_real_elem_set_error(C, 1, _mp_real_err_exp_clamp(2 * e, FLINT_BITS, prec));
        }
        else
        {
            /* expm1(m) and expm1(-m) have opposite signs */
            mp_real_init(A);
            mp_real_init(Bn);
            mp_real_init(one);
            mp_real_expm1_bits(A, &mid, prec + 4);
            mp_real_neg(Bn, &mid);
            mp_real_expm1_bits(Bn, Bn, prec + 4);
            if (rs != NULL)
            {
                mp_real_sub(S, A, Bn, wp);
                mp_real_mul_2exp_si(S, S, -1);
            }
            if (rc != NULL)
            {
                mp_real_set_ui(one, 1);
                mp_real_add(C, A, Bn, wp);
                mp_real_mul_2exp_si(C, C, -1);
                mp_real_add(C, C, one, wp);
            }
            mp_real_clear(A);
            mp_real_clear(Bn);
            mp_real_clear(one);
        }
    }

    /* the radius: over the ball, |y| <= t = |m| + rho and |sinh'| =
       cosh(y) <= exp(t), |cosh'| = |sinh(y)| <= sinh(t) <= min(t, 1)
       exp(t); for rho >= 1, the differences |sinh y - sinh m| <=
       2 sinh t and |cosh y - cosh m| <= cosh t are below exp(t) */
    if (xerr != 0)
    {
        mp_real_t R, T, E;

        mp_real_init(R);
        mp_real_init(T);
        mp_real_init(E);
        _mp_real_set_mpn_2exp(R, &xerr, 1, FLINT_BITS * xanc);
        mp_real_set(T, &mid);
        T->negative = 0;
        mp_real_add(T, T, R, 2);
        mp_real_exp_bits(E, T, 30);
        if (mp_real_abs_bound_lt_2exp_si(R) <= 0)
            mp_real_mul(E, E, R, 2);
        else
            R->size = 0;
        if (rs != NULL)
            _mp_real_elem_add_mag(S, E);
        if (rc != NULL)
        {
            if (R->size != 0 && mp_real_abs_bound_lt_2exp_si(T) <= 0)
                mp_real_mul(E, E, T, 2);
            _mp_real_elem_add_mag(C, E);
        }
        mp_real_clear(R);
        mp_real_clear(T);
        mp_real_clear(E);
    }

    if (rs != NULL)
        mp_real_swap(rs, S);
    if (rc != NULL)
        mp_real_swap(rc, C);
    mp_real_clear(S);
    mp_real_clear(C);
}

/* tanh **********************************************************************/

void
mp_real_tanh_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    mp_real_t T, u, two;
    ulong xerr;
    slong xanc, wp;

    prec = FLINT_MAX(prec, 2);
    wp = _ec_wp(prec);
    EC_SPLIT(mid, x, xerr, xanc);

    mp_real_init(T);

    if (mid.size != 0)
    {
        slong emid = _ec_emid(&mid);
        ulong L;

        /* floor |m| (capped) */
        if (mid.exp >= 2)
            L = UWORD(1) << (FLINT_BITS - 6);
        else if (mid.exp == 1)
            L = FLINT_MIN(mid.d[mid.size - 1], UWORD(1) << (FLINT_BITS - 6));
        else
            L = 0;

        if (emid <= -(prec / 2 + 4))
        {
            /* |tanh m - m| <= |m|^3 / 3 */
            mp_real_set(T, &mid);
            mp_real_add_error_2exp_si(T, _mp_real_err_exp_clamp(3 * emid, emid + FLINT_BITS, prec));
        }
        else if (L >= (ulong) prec / 2 + 4)
        {
            /* 1 - tanh |m| = 2 / (exp(2|m|) + 1) < 2^(1 - 2L) */
            _mp_real_elem_set_error(T, 1, _mp_real_err_exp_clamp(1 - 2 * (slong) L, FLINT_BITS, prec));
            if (mid.negative)
                mp_real_neg(T, T);
        }
        else
        {
            mp_real_init(u);
            mp_real_init(two);
            mp_real_set_ui(two, 2);
            mp_real_mul_2exp_si(u, &mid, 1);
            if (mid.negative)
                mp_real_neg(u, u);
            mp_real_expm1_bits(u, u, prec + 8);
            mp_real_add(T, u, two, wp);
            mp_real_div(T, u, T, wp);
            if (mid.negative)
                mp_real_neg(T, T);
            mp_real_clear(u);
            mp_real_clear(two);
        }
    }

    if (xerr != 0)
    {
        mp_real_t R, D;
        slong k = 0, kmax = WORD(1) << (FLINT_BITS - 4);

        mp_real_init(R);
        mp_real_init(D);
        _mp_real_set_mpn_2exp(R, &xerr, 1, FLINT_BITS * xanc);

        /* k <= min(2 log2(e) |y|, kmax) over the ball, from D = |m| - rho:
           in double precision (rounded down), or for D beyond the double
           range, D >= 2^(emid(D) - 2) */
        if (mid.size != 0)
        {
            mp_real_set(D, &mid);
            D->negative = 0;
            mp_real_sub(D, D, R, 2);

            if (D->size != 0 && !D->negative)
            {
                if (_ec_emid(D) > 1000)
                {
                    if (mp_real_rel_radius_lt_2exp_si(D) <= -2)
                        k = kmax;
                }
                else
                {
                    double md, mr, lo;
                    mp_real_get_dfloat(&md, &mr, 1, D);
                    lo = (md - mr) * 2.8853900817779268 * (1.0 - 0x1p-40);
                    if (lo >= (double) kmax)
                        k = kmax;
                    else if (lo >= 1.0)
                        k = (slong) lo;
                }
            }
        }

        /* for k >= 1: tanh'(y) = sech(y)^2 < 4 exp(-2|y|) <= 2^(2 - k),
           and for rho >= 1, |tanh y - tanh m| < 2 exp(-2|y|) <= 2^(1 - k)
           (y on the side of m); else |tanh'| <= 1 */
        if (k >= 1)
        {
            if (mp_real_abs_bound_lt_2exp_si(R) <= 0)
                mp_real_mul_2exp_si(R, R, 2 - k);
            else
            {
                mp_real_set_ui(R, 1);
                mp_real_mul_2exp_si(R, R, 1 - k);
            }
            /* (the radius clamped at 2^-(prec + FLINT_BITS)) */
            if (mp_real_abs_bound_lt_2exp_si(R) < -prec - FLINT_BITS)
            {
                mp_real_set_ui(R, 1);
                mp_real_mul_2exp_si(R, R, -prec - FLINT_BITS);
            }
            _mp_real_elem_add_mag(T, R);
        }
        else
            _mp_real_elem_add_rad(T, xerr, xanc);

        mp_real_clear(R);
        mp_real_clear(D);
    }

    mp_real_swap(res, T);
    mp_real_clear(T);
}

/* asin, acos ******************************************************************/

/* which = 0: asin, 1: acos */
static int
_ec_asin_acos(mp_real_t res, const mp_real_t x, slong prec, int which)
{
    mp_real_struct mid;
    mp_real_t Y, a, b, t, one;
    ulong xerr;
    slong xanc, wp, ns;
    int c, ok = 1;

    prec = FLINT_MAX(prec, 2);
    wp = _ec_wp(prec);
    EC_SPLIT(mid, x, xerr, xanc);

    c = _ec_cmpabs_one(&mid);
    if (c > 0 || (c == 0 && xerr != 0))
    {
        mp_real_zero(res);
        return 0;
    }

    mp_real_init(Y);
    mp_real_init(a);
    mp_real_init(b);
    mp_real_init(t);
    mp_real_init(one);
    mp_real_set_ui(one, 1);
    ns = _ec_exact_limbs(&mid, wp);

    if (c == 0)
    {
        /* asin(+-1) = +-pi/2, acos(1) = 0, acos(-1) = pi */
        if (which == 0 || mid.negative)
        {
            mp_real_const_pi4(Y, wp, 1);
            mp_real_mul_2exp_si(Y, Y, (which == 0) ? 1 : 2);
            if (which == 0 && mid.negative)
                mp_real_neg(Y, Y);
        }
        else
            mp_real_zero(Y);
    }
    else if (which == 0)
    {
        if (mid.size == 0)
            mp_real_zero(Y);
        else
        {
            /* atan(|m| / sqrt((1 - |m|)(1 + |m|))) */
            mp_real_set(t, &mid);
            t->negative = 0;
            mp_real_sub(a, one, t, ns);
            mp_real_add(b, one, t, ns);
            mp_real_mul(b, a, b, wp);
            mp_real_sqrt(b, b, wp);
            mp_real_div(t, t, b, wp);
            mp_real_atan_bits(Y, t, prec + 8);
            if (mid.negative)
                mp_real_neg(Y, Y);
        }
    }
    else
    {
        /* 2 atan(sqrt((1 - m) / (1 + m))) */
        mp_real_sub(a, one, &mid, ns);
        mp_real_add(b, one, &mid, ns);
        mp_real_div(t, a, b, wp);
        mp_real_sqrt(t, t, wp);
        mp_real_atan_bits(Y, t, prec + 8);
        mp_real_mul_2exp_si(Y, Y, 1);
    }

    /* |f'(y)| <= 1 / sqrt(1 - |y|) <= 1 / sqrt(1 - |m| - rho); and for a
       ball within [0, 1] (or [-1, 0]), |f(y) - f(m)| <= acos(|m| - rho)
       <= 2 sqrt(1 - |m| + rho) */
    if (xerr != 0)
    {
        mp_real_set(t, &mid);
        t->negative = 0;
        mp_real_sub(a, one, t, ns);
        ok = _ec_add_rad_div(Y, a, xerr, xanc, 1, 2, ns);
    }

    if (ok)
        mp_real_swap(res, Y);
    else
        mp_real_zero(res);

    mp_real_clear(Y);
    mp_real_clear(a);
    mp_real_clear(b);
    mp_real_clear(t);
    mp_real_clear(one);
    return ok;
}

int
mp_real_asin_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    return _ec_asin_acos(res, x, prec, 0);
}

int
mp_real_acos_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    return _ec_asin_acos(res, x, prec, 1);
}

/* asinh, acosh, atanh **********************************************************/

void
mp_real_asinh_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    mp_real_t Y, q, h, one;
    ulong xerr;
    slong xanc, wp;

    prec = FLINT_MAX(prec, 2);
    wp = _ec_wp(prec);
    EC_SPLIT(mid, x, xerr, xanc);

    mp_real_init(Y);

    if (mid.size != 0)
    {
        slong emid = _ec_emid(&mid);

        if (emid <= -(prec / 2 + 4))
        {
            /* |asinh m - m| <= |m|^3 / 6 */
            mp_real_set(Y, &mid);
            mp_real_add_error_2exp_si(Y, _mp_real_err_exp_clamp(3 * emid, emid + FLINT_BITS, prec));
        }
        else if (emid >= prec / 2 + 4)
        {
            /* asinh |m| = log(2|m|) + delta, 0 <= delta <= 1/(4 m^2)
               <= 2^(2 - 2 emid) */
            mp_real_mul_2exp_si(Y, &mid, 1);
            Y->negative = 0;
            mp_real_log_bits(Y, Y, prec + 4);
            mp_real_add_error_2exp_si(Y, _mp_real_err_exp_clamp(2 - 2 * emid, FLINT_BITS, prec));
            if (mid.negative)
                mp_real_neg(Y, Y);
        }
        else
        {
            /* log1p(|m| + m^2 / (1 + sqrt(1 + m^2))) */
            mp_real_init(q);
            mp_real_init(h);
            mp_real_init(one);
            mp_real_set_ui(one, 1);
            mp_real_mul(q, &mid, &mid, wp);
            mp_real_add(h, q, one, wp);
            mp_real_sqrt(h, h, wp);
            mp_real_add(h, h, one, wp);
            mp_real_div(q, q, h, wp);
            mp_real_set(h, &mid);
            h->negative = 0;
            mp_real_add(q, q, h, wp);
            mp_real_log1p_bits(Y, q, prec + 4);
            if (mid.negative)
                mp_real_neg(Y, Y);
            mp_real_clear(q);
            mp_real_clear(h);
            mp_real_clear(one);
        }
    }

    /* |asinh'(y)| = 1 / sqrt(1 + y^2) <= min(1, 1 / (|m| - rho)) */
    if (xerr != 0)
    {
        int done = 0;

        if (mid.size != 0 && _ec_emid(&mid) >= 2)
        {
            mp_real_struct am = mid;
            am.negative = 0;
            done = _ec_add_rad_div(Y, &am, xerr, xanc, 0, 0, wp);
        }

        if (!done)
            _mp_real_elem_add_rad(Y, xerr, xanc);
    }

    mp_real_swap(res, Y);
    mp_real_clear(Y);
}

int
mp_real_acosh_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    mp_real_t Y, d, p, one;
    ulong xerr;
    slong xanc, wp, ns;
    int c, ok = 1;

    prec = FLINT_MAX(prec, 2);
    wp = _ec_wp(prec);
    EC_SPLIT(mid, x, xerr, xanc);

    c = _ec_cmpabs_one(&mid);
    if (mid.size == 0 || mid.negative || c < 0 || (c == 0 && xerr != 0))
    {
        mp_real_zero(res);
        return 0;
    }

    mp_real_init(Y);
    mp_real_init(d);
    mp_real_init(p);
    mp_real_init(one);
    mp_real_set_ui(one, 1);
    ns = _ec_exact_limbs(&mid, wp);

    /* d = m - 1, exact */
    mp_real_sub(d, &mid, one, ns);

    if (c == 0)
        mp_real_zero(Y);
    else
    {
        slong emid = _ec_emid(&mid);

        if (emid >= prec / 2 + 4)
        {
            /* acosh m = log(2m) - delta, 0 <= delta <= 1/m^2 <= 2^(2 - 2 emid) */
            mp_real_mul_2exp_si(Y, &mid, 1);
            mp_real_log_bits(Y, Y, prec + 4);
            mp_real_add_error_2exp_si(Y, _mp_real_err_exp_clamp(2 - 2 * emid, FLINT_BITS, prec));
        }
        else
        {
            /* log1p(d + sqrt(d (m + 1))) */
            mp_real_add(p, &mid, one, wp);
            mp_real_mul(p, p, d, wp);
            mp_real_sqrt(p, p, wp);
            mp_real_add(p, p, d, wp);
            mp_real_log1p_bits(Y, p, prec + 4);
        }
    }

    /* |acosh'(y)| = 1 / sqrt((y - 1)(y + 1)) <= 1 / sqrt(D (D + 2)) with
       D = m - 1 - rho, which is at most 1 / sqrt(D) and at most 1 / D
       (the better one for d = m - 1 >= 2); and |acosh y - acosh m| <=
       acosh(1 + d + rho) <= sqrt(2 (d + rho)) */
    if (xerr != 0)
        ok = _ec_add_rad_div(Y, d, xerr, xanc, (d->size == 0 || _ec_emid(d) <= 1), 1, ns);

    if (ok)
        mp_real_swap(res, Y);
    else
        mp_real_zero(res);

    mp_real_clear(Y);
    mp_real_clear(d);
    mp_real_clear(p);
    mp_real_clear(one);
    return ok;
}

int
mp_real_atanh_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    mp_real_t Y, a, t, one;
    ulong xerr;
    slong xanc, wp, ns;
    int ok = 1;

    prec = FLINT_MAX(prec, 2);
    wp = _ec_wp(prec);
    EC_SPLIT(mid, x, xerr, xanc);

    if (_ec_cmpabs_one(&mid) >= 0)
    {
        mp_real_zero(res);
        return 0;
    }

    mp_real_init(Y);
    mp_real_init(a);
    mp_real_init(t);
    mp_real_init(one);
    mp_real_set_ui(one, 1);
    ns = _ec_exact_limbs(&mid, wp);

    /* a = 1 - |m|, exact */
    mp_real_set(t, &mid);
    t->negative = 0;
    mp_real_sub(a, one, t, ns);

    if (mid.size != 0)
    {
        /* log1p(2|m| / (1 - |m|)) / 2 */
        mp_real_div(t, t, a, wp);
        mp_real_mul_2exp_si(t, t, 1);
        mp_real_log1p_bits(Y, t, prec + 4);
        mp_real_mul_2exp_si(Y, Y, -1);
        if (mid.negative)
            mp_real_neg(Y, Y);
    }

    /* |atanh'(y)| = 1 / (1 - y^2) <= 1 / (1 - |m| - rho) */
    if (xerr != 0)
        ok = _ec_add_rad_div(Y, a, xerr, xanc, 0, 0, ns);

    if (ok)
        mp_real_swap(res, Y);
    else
        mp_real_zero(res);

    mp_real_clear(Y);
    mp_real_clear(a);
    mp_real_clear(t);
    mp_real_clear(one);
    return ok;
}
