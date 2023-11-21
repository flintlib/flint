/*
    Copyright (C) 2014, 2018, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* We do Newton iteration in floating-point arithmetic with some
   guard bits. With GUARD_BITS = 32 the result is certainly accurate
   to more than GUARD_BITS / 2 extra bits. With a detailed error
   analysis this could be set much more conservatively
   (see e.g. MPFR's algorithms documentation), but currently
   we only care about extremely high precision and we always
   start from an extremely high precision initial value,
   so fine-tuning the precision management isn't going to help
   performance a whole lot. */
#define GUARD_BITS 32
#define INV_NEWTON_CUTOFF 24000
#define DIV_NEWTON_CUTOFF 70000

#ifdef FLINT_HAVE_FFT_SMALL
#define WANT_NEWTON(prec, xbits, ybits) (prec) >= INV_NEWTON_CUTOFF && (ybits) > (prec) * 0.5 && ((xbits) < (prec) * 0.01 || (prec) >= DIV_NEWTON_CUTOFF)
#else
#define WANT_NEWTON(prec, xbits, ybits) 0
#endif

void
_arf_inv_newton(arf_t res, const arf_t x, slong prec)
{
    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    if (prec < INV_NEWTON_CUTOFF)
    {
        arf_set_round(res, x, wp, ARF_RND_DOWN);
        arf_ui_div(res, 1, res, wp, ARF_RND_DOWN);
    }
    else
    {
        arf_t r, t;

        arf_init(r);
        arf_init(t);

        _arf_inv_newton(r, x, hp);

        /* r - r*(x*r - 1) */

        if (arf_bits(x) <= wp)
        {
            arf_mul(t, x, r, wp, ARF_RND_DOWN);
        }
        else
        {
            arf_set_round(t, x, wp, ARF_RND_DOWN);
            arf_mul(t, t, r, wp, ARF_RND_DOWN);
        }

        arf_sub_ui(t, t, 1, hp, ARF_RND_DOWN);
        arf_mul(t, t, r, hp, ARF_RND_DOWN);
        arf_sub(res, r, t, wp, ARF_RND_DOWN);

        arf_clear(r);
        arf_clear(t);
    }
}

/* Karp-Markstein */
void
_arf_div_newton(arf_t res, const arf_t x, const arf_t y, slong prec)
{
    arf_t xn, yn, t;

    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    arf_init(xn);
    arf_init(yn);
    arf_init(t);

    _arf_inv_newton(xn, y, hp);
    arf_set_round(t, x, hp, ARF_RND_DOWN);
    arf_mul(yn, xn, t, hp, ARF_RND_DOWN);
    arf_mul(t, y, yn, wp, ARF_RND_DOWN);
    arf_sub(t, x, t, hp, ARF_RND_DOWN);
    arf_mul(t, t, xn, hp, ARF_RND_DOWN);
    arf_add(res, yn, t, wp, ARF_RND_DOWN);

    arf_clear(xn);
    arf_clear(yn);
    arf_clear(t);
}

void
arb_div_newton(arb_t res, const arb_t x, const arb_t y, slong prec)
{
    mag_t zr, xm, ym, yl, yw;

    if (arf_is_special(arb_midref(x)) || arf_is_special(arb_midref(y)))
    {
        arb_indeterminate(res);
        return;
    }

    mag_init_set_arf(xm, arb_midref(x));
    mag_init_set_arf(ym, arb_midref(y));
    mag_init(zr);
    mag_init(yl);
    mag_init(yw);

    /* (|x|*yrad + |y|*xrad)/(|y|*(|y|-yrad)) */
    mag_mul(zr, xm, arb_radref(y));
    mag_addmul(zr, ym, arb_radref(x));
    arb_get_mag_lower(yw, y);
    arf_get_mag_lower(yl, arb_midref(y));
    mag_mul_lower(yl, yl, yw);
    mag_div(zr, zr, yl);

    _arf_div_newton(arb_midref(res), arb_midref(x), arb_midref(y), prec);
    arf_mag_add_ulp(arb_radref(res), zr, arb_midref(res), prec + GUARD_BITS / 2);
    arb_set_round(res, res, prec);

    mag_clear(xm);
    mag_clear(ym);
    mag_clear(zr);
    mag_clear(yl);
    mag_clear(yw);
}

void
arb_div_arf_newton(arb_t res, const arb_t x, const arf_t y, slong prec)
{
    mag_t zr, ym;

    if (arf_is_special(arb_midref(x)) || arf_is_special(y))
    {
        arb_indeterminate(res);
        return;
    }

    mag_init(ym);
    mag_init(zr);

    arf_get_mag_lower(ym, y);
    mag_div(zr, arb_radref(x), ym);

    _arf_div_newton(arb_midref(res), arb_midref(x), y, prec);
    arf_mag_add_ulp(arb_radref(res), zr, arb_midref(res), prec + GUARD_BITS / 2);
    arb_set_round(res, res, prec);

    mag_clear(ym);
    mag_clear(zr);
}

void
arb_div_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    mag_t zr, ym;
    int inexact;

    if (arf_is_special(y) || !arb_is_finite(x))
    {
        if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)) &&
            arf_is_finite(y) && !arf_is_zero(y))
        {
            /* inf / finite nonzero = inf */
            arf_div(arb_midref(z), arb_midref(x), y, prec, ARF_RND_DOWN);
            mag_zero(arb_radref(z));
        }
        else if (arb_is_finite(x) && arf_is_inf(y))
        {
            /* finite / inf = 0 */
            arb_zero(z);
        }
        else if (!arf_is_nan(arb_midref(x)) && mag_is_inf(arb_radref(x)) &&
            arf_is_finite(y) && !arf_is_zero(y))
        {
            /* [+/- inf] / finite nonzero = [+/- inf] */
            arb_zero_pm_inf(z);
        }
        else
        {
            arb_indeterminate(z);
        }
    }
    else if (mag_is_zero(arb_radref(x)))
    {
        inexact = arf_div(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
    else if (WANT_NEWTON(prec, arb_bits(x), arf_bits(y)))
    {
        arb_div_arf_newton(z, x, y, prec);
    }
    else
    {
        mag_init(ym);
        mag_init(zr);

        arf_get_mag_lower(ym, y);
        mag_div(zr, arb_radref(x), ym);

        inexact = arf_div(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(ym);
        mag_clear(zr);
    }
}

static void
arb_div_wide(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    mag_t a, b, t, u;

    mag_init(t);
    arb_get_mag_lower(t, y);

    if (mag_is_zero(t))
    {
        arb_indeterminate(z);
    }
    else if (arf_is_zero(arb_midref(x)))
    {
        mag_div(arb_radref(z), arb_radref(x), t);
        arf_zero(arb_midref(z));
    }
    else
    {
        if (arf_cmpabs_mag(arb_midref(x), arb_radref(x)) >= 0)
        {
            /* [a,b] /  [t,u] =  [a/u, b/t]
               [a,b] / -[t,u] = -[a/u, b/t]
              -[a,b] /  [t,u] = -[a/u, b/t]
              -[a,b] / -[t,u] =  [a/u, b/t] */
            mag_init(a);
            mag_init(b);
            mag_init(u);

            arb_get_mag_lower(a, x);
            arb_get_mag(b, x);
            arb_get_mag(u, y);

            mag_div_lower(a, a, u);
            mag_div(b, b, t);

            if ((arf_sgn(arb_midref(x)) < 0) ^ (arf_sgn(arb_midref(y)) < 0))
            {
                arb_set_interval_mag(z, a, b, prec);
                arb_neg(z, z);
            }
            else
            {
                arb_set_interval_mag(z, a, b, prec);
            }

            mag_clear(a);
            mag_clear(b);
            mag_clear(u);
        }
        else
        {
            /* [-a,b] /  [t,u] = [-a/t, b/t]
               [-a,b] / -[t,u] = [-b/t, a/t] */
            mag_init(a);
            mag_init(b);

            arb_get_mag(b, x);
            arf_get_mag_lower(a, arb_midref(x));
            mag_sub(a, arb_radref(x), a);

            if ((arf_sgn(arb_midref(x)) < 0) ^ (arf_sgn(arb_midref(y)) < 0))
                mag_swap(a, b);

            mag_div(a, a, t);
            mag_div(b, b, t);
            arb_set_interval_neg_pos_mag(z, a, b, prec);

            mag_clear(a);
            mag_clear(b);
        }
    }

    mag_clear(t);
}

void
arb_div(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    int inexact;
    slong acc, xacc, yacc;

    if (mag_is_zero(arb_radref(y)))
    {
        arb_div_arf(z, x, arb_midref(y), prec);
    }
    else if (arf_is_zero(arb_midref(y))) /* anything / 0 = nan */
    {
        arb_indeterminate(z);
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y)) /* fast case */
    {
        /* both finite, both midpoint and radius of y not special */
        yacc = MAG_EXP(arb_midref(y)) - ARF_EXP(arb_radref(y));
        xacc = arb_rel_accuracy_bits(x);
        acc = FLINT_MIN(xacc, yacc);
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc <= 20)
        {
            arb_div_wide(z, x, y, prec);
        }
        else if (WANT_NEWTON(prec, arb_bits(x), arb_bits(y)))
        {
            arb_div_newton(z, x, y, prec);
        }
        else
        {
            mag_t t, u, v;
            mag_init(t);  /* no need to free */
            mag_init(u);
            mag_init(v);

            /*               (x*yrad + y*xrad)/(y*(y-yrad))
                 <=  (1+eps) (x*yrad + y*xrad)/y^2
                 <=  (1+eps) (x*yrad/y^2 + y*xrad/y^2)
                 <=  (1+eps) ((x/y)*yrad/y + xrad/y)
                 <=  (1+eps) ((x/y)*yrad + xrad)/y
            */

            /* t = y */
            arf_get_mag_lower(t, arb_midref(y));

            inexact = arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

            /* u = x/y */
            arf_get_mag(u, arb_midref(z));

            /* v = (x/y)*yrad + xrad */
            MAG_MAN(v) = MAG_MAN(arb_radref(x));
            MAG_EXP(v) = MAG_EXP(arb_radref(x));
            mag_fast_addmul(v, arb_radref(y), u);
            /* v = v / y */
            mag_div(arb_radref(z), v, t);

            /* correct for errors */
            MAG_MAN(t) = MAG_ONE_HALF + (MAG_ONE_HALF >> 16);
            MAG_EXP(t) = 1;
            mag_fast_mul(arb_radref(z), arb_radref(z), t);

            if (inexact)
                arf_mag_fast_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        }
    }
    else if (!arb_is_finite(y))
    {
        /* finite / inf = 0 */
        if (arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)) &&
            arb_is_finite(x))
            arb_zero(z);
        else
            arb_indeterminate(z);
    }
    else if (!arb_is_finite(x))
    {
        if (arb_contains_zero(y) || arf_is_nan(arb_midref(x)))
        {
            arb_indeterminate(z);
        }
        else
        {
            /* +/- inf  /  finite nonzero  = +/- inf */
            if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            {
                arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
                mag_zero(arb_radref(z));
            }
            else if (!arf_is_nan(arb_midref(x)) && mag_is_inf(arb_radref(x)))
            {
                arb_zero_pm_inf(z);
            }
            else
            {
                arb_indeterminate(z);
            }
        }
    }
    else
    {
        yacc = arb_rel_accuracy_bits(y);
        xacc = arb_rel_accuracy_bits(x);
        acc = FLINT_MIN(xacc, yacc);
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc <= 20)
        {
            arb_div_wide(z, x, y, prec);
        }
        else
        {
            mag_t zr, xm, ym, yl, yw;

            mag_init_set_arf(xm, arb_midref(x));
            mag_init_set_arf(ym, arb_midref(y));
            mag_init(zr);
            mag_init(yl);
            mag_init(yw);

            /* (|x|*yrad + |y|*xrad)/(|y|*(|y|-yrad)) */
            mag_mul(zr, xm, arb_radref(y));
            mag_addmul(zr, ym, arb_radref(x));
            arb_get_mag_lower(yw, y);

            arf_get_mag_lower(yl, arb_midref(y));
            mag_mul_lower(yl, yl, yw);

            mag_div(zr, zr, yl);

            inexact = arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

            if (inexact)
                arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
            else
                mag_swap(arb_radref(z), zr);

            mag_clear(xm);
            mag_clear(ym);
            mag_clear(zr);
            mag_clear(yl);
            mag_clear(yw);
        }
    }
}

void
arb_div_si(arb_t z, const arb_t x, slong y, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_ui(arb_t z, const arb_t x, ulong y, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_div_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_div_arf(z, x, t, prec);
        arf_clear(t);
    }
}

void
arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, slong prec)
{
    if (WANT_NEWTON(prec, fmpz_bits(x), fmpz_bits(y)))
    {
        arb_t t;
        arb_init(t);
        arb_set_round_fmpz(t, y, prec + 10);
        arb_set_round_fmpz(z, x, prec + 10);
        arb_div_newton(z, z, t, prec);
        arb_clear(t);
    }
    else
    {
        int inexact = arf_fmpz_div_fmpz(arb_midref(z), x, y, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
}

void
arb_ui_div(arb_t z, ulong x, const arb_t y, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, x);
    arb_div(z, t, y, prec);
    arb_clear(t);
}

/* important: must be a valid divexact */
void
_arb_fmpz_divapprox_newton(fmpz_t res, const fmpz_t x, const fmpz_t y, slong exp)
{
    slong xb, yb, zb, prec;
    arf_t t, u;

    xb = fmpz_bits(x);
    yb = fmpz_bits(y);
    zb = xb - yb + exp;
    prec = FLINT_MAX(zb, 0) + 16;

    arf_init(t);
    arf_init(u);

    arf_set_round_fmpz(t, x, prec, ARF_RND_NEAR);
    arf_mul_2exp_si(t, t, exp);
    arf_set_round_fmpz(u, y, prec, ARF_RND_NEAR);
    _arf_div_newton(t, t, u, prec);
    arf_get_fmpz(res, t, ARF_RND_NEAR);

    arf_clear(t);
    arf_clear(u);
}

void
arb_fmpz_divapprox(fmpz_t res, const fmpz_t x, const fmpz_t y)
{
#ifdef FLINT_HAVE_FFT_SMALL
    if (!COEFF_IS_MPZ(*x) || !COEFF_IS_MPZ(*y))
    {
        fmpz_tdiv_q(res, x, y);
    }
    else
    {
        slong xb, yb;

        xb = mpz_size(COEFF_TO_PTR(*x));
        yb = mpz_size(COEFF_TO_PTR(*y));

        xb *= FLINT_BITS;
        yb *= FLINT_BITS;

        if (xb - yb >= DIV_NEWTON_CUTOFF && yb >= DIV_NEWTON_CUTOFF)
        {
            _arb_fmpz_divapprox_newton(res, x, y, 0);
        }
        else
        {
            fmpz_tdiv_q(res, x, y);
        }
    }
#else
    fmpz_tdiv_q(res, x, y);
#endif
}
