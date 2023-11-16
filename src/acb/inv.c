/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

static void
_arb_arf_div_rounded_den(arb_t res, const arf_t x, const arf_t y, int yinexact, slong prec)
{
    int inexact = arf_div(arb_midref(res), x, y, prec, ARB_RND);

    if (yinexact && !arf_is_special(arb_midref(res)))
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec - 1);
    else if (inexact)
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
    else
        mag_zero(arb_radref(res));
}

static void
_arb_arf_div_rounded_den_add_err(arb_t res, const arf_t x, const arf_t y, int yinexact, slong prec)
{
    int inexact = arf_div(arb_midref(res), x, y, prec, ARB_RND);

    if (yinexact && !arf_is_special(arb_midref(res)))
        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec - 1);
    else if (inexact)
        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);
}

void
acb_inv(acb_t res, const acb_t z, slong prec)
{
    mag_t am, bm;
    slong hprec;

#define a arb_midref(acb_realref(z))
#define b arb_midref(acb_imagref(z))
#define x arb_radref(acb_realref(z))
#define y arb_radref(acb_imagref(z))

    /* choose precision for the floating-point approximation of a^2+b^2 so
       that the double rounding result in less than
       2 ulp error; also use at least MAG_BITS bits since the
       value will be recycled for error bounds */
    hprec = FLINT_MAX(prec + 3, MAG_BITS);

    if (arb_is_zero(acb_imagref(z)))
    {
        arb_inv(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
        return;
    }

    if (arb_is_zero(acb_realref(z)))
    {
        arb_inv(acb_imagref(res), acb_imagref(z), prec);
        arb_neg(acb_imagref(res), acb_imagref(res));
        arb_zero(acb_realref(res));
        return;
    }

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (mag_is_zero(x) && mag_is_zero(y))
    {
        int inexact;

        arf_t a2b2;
        arf_init(a2b2);

        inexact = arf_sosq(a2b2, a, b, hprec, ARF_RND_DOWN);

        if (arf_is_special(a2b2))
        {
            acb_indeterminate(res);
        }
        else
        {
            _arb_arf_div_rounded_den(acb_realref(res), a, a2b2, inexact, prec);
            _arb_arf_div_rounded_den(acb_imagref(res), b, a2b2, inexact, prec);
            arf_neg(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(res)));
        }

        arf_clear(a2b2);
        return;
    }

    mag_init(am);
    mag_init(bm);

    /* first bound |a|-x, |b|-y */
    arb_get_mag_lower(am, acb_realref(z));
    arb_get_mag_lower(bm, acb_imagref(z));

    if ((mag_is_zero(am) && mag_is_zero(bm)))
    {
        acb_indeterminate(res);
    }
    else
    {
        /*
        The propagated error in the real part is given exactly by

             (a+x')/((a+x')^2+(b+y'))^2 - a/(a^2+b^2) = P / Q,

             P = [(b^2-a^2) x' - a (x'^2+y'^2 + 2y'b)]
             Q = [(a^2+b^2)((a+x')^2+(b+y')^2)]

        where |x'| <= x and |y'| <= y, and analogously for the imaginary part.
        */
        mag_t t, u, v, w;
        arf_t a2b2;
        int inexact;

        mag_init(t);
        mag_init(u);
        mag_init(v);
        mag_init(w);

        arf_init(a2b2);

        inexact = arf_sosq(a2b2, a, b, hprec, ARF_RND_DOWN);

        /* compute denominator */
        /* t = (|a|-x)^2 + (|b|-x)^2 (lower bound) */
        mag_mul_lower(t, am, am);
        mag_mul_lower(u, bm, bm);
        mag_add_lower(t, t, u);
        /* u = a^2 + b^2 (lower bound) */
        arf_get_mag_lower(u, a2b2);
        /* t = ((|a|-x)^2 + (|b|-x)^2)(a^2 + b^2) (lower bound) */
        mag_mul_lower(t, t, u);

        /* compute numerator */
        /* real: |a^2-b^2| x  + |a| ((x^2 + y^2) + 2 |b| y)) */
        /* imag: |a^2-b^2| y  + |b| ((x^2 + y^2) + 2 |a| x)) */
        /* am, bm = upper bounds for a, b */
        arf_get_mag(am, a);
        arf_get_mag(bm, b);

        /* v = x^2 + y^2 */
        mag_mul(v, x, x);
        mag_addmul(v, y, y);

        /* u = |a| ((x^2 + y^2) + 2 |b| y) */
        mag_mul_2exp_si(u, bm, 1);
        mag_mul(u, u, y);
        mag_add(u, u, v);
        mag_mul(u, u, am);

        /* v = |b| ((x^2 + y^2) + 2 |a| x) */
        mag_mul_2exp_si(w, am, 1);
        mag_addmul(v, w, x);
        mag_mul(v, v, bm);

        /* w = |b^2 - a^2| (upper bound) */
        if (arf_cmpabs(a, b) >= 0)
            mag_mul(w, am, am);
        else
            mag_mul(w, bm, bm);

        mag_addmul(u, w, x);
        mag_addmul(v, w, y);

        mag_div(arb_radref(acb_realref(res)), u, t);
        mag_div(arb_radref(acb_imagref(res)), v, t);

        _arb_arf_div_rounded_den_add_err(acb_realref(res), a, a2b2, inexact, prec);
        _arb_arf_div_rounded_den_add_err(acb_imagref(res), b, a2b2, inexact, prec);
        arf_neg(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(res)));

        mag_clear(t);
        mag_clear(u);
        mag_clear(v);
        mag_clear(w);

        arf_clear(a2b2);
    }

    mag_clear(am);
    mag_clear(bm);
#undef a
#undef b
#undef x
#undef y
}

