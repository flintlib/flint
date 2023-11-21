/*
    Copyright (C) 2014, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* see comments in div.c */
#define GUARD_BITS 32
#define RSQRT_NEWTON_CUTOFF 4000
#define SQRT_NEWTON_CUTOFF 200000

void
_arf_rsqrt_newton(arf_t res, const arf_t x, slong prec)
{
    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    if (prec < RSQRT_NEWTON_CUTOFF)
    {
        arf_set_round(res, x, wp, ARF_RND_DOWN);
        arf_rsqrt(res, res, wp, ARF_RND_DOWN);
    }
    else
    {
        arf_t r, t, u;

        arf_init(r);
        arf_init(t);
        arf_init(u);

        _arf_rsqrt_newton(r, x, hp);

        /* r - r*(x*r^2 - 1)/2 */

        arf_mul(t, r, r, wp, ARF_RND_DOWN);

        if (arf_bits(x) <= wp)
        {
            arf_mul(t, t, x, wp, ARF_RND_DOWN);
        }
        else
        {
            arf_set_round(u, x, wp, ARF_RND_DOWN);
            arf_mul(t, t, u, wp, ARF_RND_DOWN);
        }

        arf_sub_ui(t, t, 1, hp, ARF_RND_DOWN);
        arf_mul_2exp_si(t, t, -1);
        arf_mul(t, t, r, hp, ARF_RND_DOWN);

        arf_sub(res, r, t, wp, ARF_RND_DOWN);

        arf_clear(r);
        arf_clear(t);
        arf_clear(u);
    }
}

void
_arf_sqrt_newton(arf_t res, const arf_t x, slong prec)
{
    arf_t t, u, v;

    slong wp = prec + GUARD_BITS;
    slong hp = prec / 2 + GUARD_BITS;

    arf_init(t);
    arf_init(u);
    arf_init(v);

    _arf_rsqrt_newton(t, x, hp);

    if (arf_bits(x) <= hp)
    {
        arf_mul(v, t, x, hp, ARF_RND_DOWN);
    }
    else
    {
        arf_set_round(u, x, hp, ARF_RND_DOWN);
        arf_mul(v, t, u, hp, ARF_RND_DOWN);
    }

    arf_mul(u, v, v, wp, ARF_RND_DOWN);
    arf_sub(u, x, u, hp, ARF_RND_DOWN);
    arf_mul(u, u, t, wp, ARF_RND_DOWN);
    arf_mul_2exp_si(u, u, -1);
    arf_add(res, v, u, wp, ARF_RND_DOWN);

    arf_clear(t);
    arf_clear(u);
    arf_clear(v);
}

void
arb_rsqrt_arf_newton(arb_t res, const arf_t x, slong prec)
{
    if (arf_is_special(x) || arf_sgn(x) < 0)
    {
        arb_indeterminate(res);
        return;
    }

    /* special case: handle 2^(2n) exactly */
    if (ARF_IS_POW2(x) && fmpz_is_odd(ARF_EXPREF(x)))
    {
        arf_rsqrt(arb_midref(res), x, prec, ARF_RND_DOWN);
        mag_zero(arb_radref(res));
        return;
    }

    _arf_rsqrt_newton(arb_midref(res), x, prec);
    arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec + GUARD_BITS / 2);
    arb_set_round(res, res, prec);
}

void
arb_rsqrt_arf(arb_t res, const arf_t x, slong prec)
{
    if (arf_is_special(x) || arf_sgn(x) < 0)
    {
        if (arf_is_pos_inf(x))
            arb_zero(res);
        else
            arb_indeterminate(res);
        return;
    }

#ifdef FLINT_HAVE_FFT_SMALL
    if (prec > RSQRT_NEWTON_CUTOFF)
    {
        arb_rsqrt_arf_newton(res, x, prec);
    }
    else
#endif
    {
        int inexact;

        inexact = arf_rsqrt(arb_midref(res), x, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
        else
            mag_zero(arb_radref(res));
    }
}

void
arb_rsqrt(arb_t z, const arb_t x, slong prec)
{
    if (!arb_is_finite(x) || arf_sgn(arb_midref(x)) <= 0)
    {
        if (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_zero(z);
        else
            arb_indeterminate(z);
    }
    else if (mag_is_zero(arb_radref(x)))
    {
        arb_rsqrt_arf(z, arb_midref(x), prec);
    }
    else
    {
        mag_t t, u;
        slong acc;

        mag_init(t);

        arb_get_mag_lower(t, x);

        acc = arb_rel_accuracy_bits(x);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc <= 20)
        {
            if (mag_is_zero(t))
            {
                arb_indeterminate(z);
            }
            else
            {
                mag_init(u);
                arb_get_mag(u, x);

                mag_rsqrt(t, t);
                mag_rsqrt_lower(u, u);

                arb_set_interval_mag(z, u, t, prec);

                mag_clear(u);
            }
        }
        else
        {
            /* error bound: (1/2) (x-r)^(-3/2) * r */
            mag_init(u);

            mag_rsqrt(u, t);
            mag_div(t, u, t);
            mag_mul(t, t, arb_radref(x));
            mag_mul_2exp_si(t, t, -1);

            arb_rsqrt_arf(z, arb_midref(x), prec);
            mag_add(arb_radref(z), arb_radref(z), t);

            mag_clear(u);
        }

        mag_clear(t);
    }
}

void
arb_sqrt_arf_newton(arb_t res, const arf_t x, slong prec)
{
    if (arf_is_special(x) || arf_sgn(x) < 0)
    {
        if (arf_is_zero(x) || arf_is_pos_inf(x))
            arb_set_arf(res, x);
        else
            arb_indeterminate(res);
        return;
    }

    /* special case: handle 2^(2n) exactly */
    /* todo: detect other simple squares? */
    if (ARF_IS_POW2(x) && fmpz_is_odd(ARF_EXPREF(x)))
    {
        arf_sqrt(arb_midref(res), x, prec, ARF_RND_DOWN);
        mag_zero(arb_radref(res));
        return;
    }

    _arf_sqrt_newton(arb_midref(res), x, prec);
    arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec + GUARD_BITS / 2);
    arb_set_round(res, res, prec);
}

void
arb_sqrt_arf(arb_t res, const arf_t x, slong prec)
{
    if (arf_is_special(x) || arf_sgn(x) < 0)
    {
        if (arf_is_zero(x) || arf_is_pos_inf(x))
            arb_set_arf(res, x);
        else
            arb_indeterminate(res);
        return;
    }

#ifdef FLINT_HAVE_FFT_SMALL
    if (prec > SQRT_NEWTON_CUTOFF)
    {
        arb_sqrt_arf_newton(res, x, prec);
    }
    else
#endif
    {
        if (arf_sqrt(arb_midref(res), x, prec, ARB_RND))
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
        else
            mag_zero(arb_radref(res));
    }
}

void
arb_sqrt_newton(arb_t z, const arb_t x, slong prec)
{
    mag_t zr, rx;

    mag_init(zr);
    mag_init(rx);

    /* rx = upper bound for r / x */
    arf_get_mag_lower(rx, arb_midref(x));
    mag_div(rx, arb_radref(x), rx);

    arb_sqrt_arf_newton(z, arb_midref(x), prec);
    /* zr = upper bound for sqrt(x) */
    arb_get_mag(zr, z);

    /* propagated error:   sqrt(x) - sqrt(x-r)
                         = sqrt(x) * [1 - sqrt(1 - r/x)]
                        <= sqrt(x) * 0.5 * (rx + rx^2)  */
    mag_addmul(rx, rx, rx);
    mag_mul(zr, zr, rx);
    mag_mul_2exp_si(zr, zr, -1);

    mag_add(arb_radref(z), arb_radref(z), zr);

    mag_clear(zr);
    mag_clear(rx);
}

void
arb_sqrt(arb_t z, const arb_t x, slong prec)
{
    mag_t rx, zr;
    int inexact;

    if (mag_is_zero(arb_radref(x)))
    {
        arb_sqrt_arf(z, arb_midref(x), prec);
    }
    else if (arf_is_special(arb_midref(x)) ||
              arf_sgn(arb_midref(x)) < 0 || mag_is_inf(arb_radref(x)))
    {
        if (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_sqrt_arf(z, arb_midref(x), prec);
        else
            arb_indeterminate(z);
    }
    else  /* now both mid and rad are non-special values, mid > 0 */
    {
        slong acc;

        acc = _fmpz_sub_small(ARF_EXPREF(arb_midref(x)), MAG_EXPREF(arb_radref(x)));
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc < 0)
        {
            arb_indeterminate(z);
        }
        else if (acc <= 20)
        {
            mag_t t, u;

            mag_init(t);
            mag_init(u);

            arb_get_mag_lower(t, x);

            if (mag_is_zero(t) && arb_contains_negative(x))
            {
                arb_indeterminate(z);
            }
            else
            {
                arb_get_mag(u, x);
                mag_sqrt_lower(t, t);
                mag_sqrt(u, u);
                arb_set_interval_mag(z, t, u, prec);
            }

            mag_clear(t);
            mag_clear(u);
        }
#ifdef FLINT_HAVE_FFT_SMALL
        else if (prec > SQRT_NEWTON_CUTOFF)
        {
            arb_sqrt_newton(z, x, prec);
        }
#endif
        else if (ARB_IS_LAGOM(x)) /* small exponents, acc *and* prec >= 20 */
        {
            mag_t t;
            mag_init(t); /* no need to free */

            inexact = arf_sqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

            /* sqrt(x) - sqrt(x-r) <= 0.5 * r * rsqrt(x-r)  */
            /* we have rsqrt(x-r) ~= 1/sqrt(x) */
            arf_get_mag_lower(t, arb_midref(z));

            /* note: we need to write rad(z) first to use fast_mul later */
            mag_div(arb_radref(z), arb_radref(x), t);

            /* We are guaranteed to have acc and prec >= 20. */
            /* 0.5 + eps corrects for errors */
            MAG_MAN(t) = MAG_ONE_HALF + (MAG_ONE_HALF >> 16);
            MAG_EXP(t) = 0;
            mag_fast_mul(arb_radref(z), arb_radref(z), t);

            if (inexact)
                arf_mag_fast_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        }
        else
        {
            mag_init(zr);
            mag_init(rx);

            /* rx = upper bound for r / x */
            arf_get_mag_lower(rx, arb_midref(x));
            mag_div(rx, arb_radref(x), rx);

            inexact = arf_sqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

            /* zr = upper bound for sqrt(x) */
            arf_get_mag(zr, arb_midref(z));
            if (inexact)
                arf_mag_add_ulp(zr, zr, arb_midref(z), prec);

            /* propagated error:   sqrt(x) - sqrt(x-r)
                                 = sqrt(x) * [1 - sqrt(1 - r/x)]
                                <= sqrt(x) * 0.5 * (rx + rx^2)  */
            mag_addmul(rx, rx, rx);
            mag_mul(zr, zr, rx);
            mag_mul_2exp_si(zr, zr, -1);

            /* add the rounding error */
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
            else
                mag_swap(arb_radref(z), zr);

            mag_clear(zr);
            mag_clear(rx);
        }
    }
}

void
arb_rsqrt_ui(arb_t z, ulong x, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    arb_rsqrt_arf(z, t, prec);
}

void
arb_sqrt_ui(arb_t z, ulong x, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    arb_sqrt_arf(z, t, prec);
}

void
arb_sqrt_fmpz(arb_t z, const fmpz_t x, slong prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    arb_sqrt_arf(z, t, prec);
    arf_clear(t);
}
