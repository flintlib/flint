/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

/* we don't need any error bounding, so we define a few helper
   functions that ignore the radii */

static inline void
acb_sub_mid(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    arf_sub(arb_midref(acb_realref(z)),
        arb_midref(acb_realref(x)),
        arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_sub(arb_midref(acb_imagref(z)),
        arb_midref(acb_imagref(x)),
        arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static inline void
acb_add_mid(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    arf_add(arb_midref(acb_realref(z)),
        arb_midref(acb_realref(x)),
        arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_add(arb_midref(acb_imagref(z)),
        arb_midref(acb_imagref(x)),
        arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static inline void
acb_mul_mid(acb_t z, const acb_t x, const acb_t y, slong prec)
{
#define a arb_midref(acb_realref(x))
#define b arb_midref(acb_imagref(x))
#define c arb_midref(acb_realref(y))
#define d arb_midref(acb_imagref(y))
#define e arb_midref(acb_realref(z))
#define f arb_midref(acb_imagref(z))

    arf_complex_mul(e, f, a, b, c, d, prec, ARF_RND_DOWN);

#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
}

static inline void
acb_inv_mid(acb_t z, const acb_t x, slong prec)
{
    arf_t t;
    arf_init(t);

#define a arb_midref(acb_realref(x))
#define b arb_midref(acb_imagref(x))
#define e arb_midref(acb_realref(z))
#define f arb_midref(acb_imagref(z))

    arf_mul(t, a, a, prec, ARF_RND_DOWN);
    arf_addmul(t, b, b, prec, ARF_RND_DOWN);

    arf_div(e, a, t, prec, ARF_RND_DOWN);
    arf_div(f, b, t, prec, ARF_RND_DOWN);

    arf_neg(f, f);

#undef a
#undef b
#undef e
#undef f

    arf_clear(t);
}

void
_acb_poly_evaluate_mid(acb_t res, acb_srcptr f, slong len,
                           const acb_t a, slong prec)
{
    slong i = len - 1;
    acb_t t;

    acb_init(t);
    acb_set(res, f + i);

    for (i = len - 2; i >= 0; i--)
    {
        acb_mul_mid(t, res, a, prec);
        acb_add_mid(res, f + i, t, prec);
    }

    acb_clear(t);
}

void
_acb_poly_refine_roots_durand_kerner(acb_ptr roots,
        acb_srcptr poly, slong len, slong prec)
{
    slong i, j;

    acb_t x, y, t;

    acb_init(x);
    acb_init(y);
    acb_init(t);

    for (i = 0; i < len - 1; i++)
    {
        _acb_poly_evaluate_mid(x, poly, len, roots + i, prec);

        acb_set(y, poly + len - 1);

        for (j = 0; j < len - 1; j++)
        {
            if (i != j)
            {
                acb_sub_mid(t, roots + i, roots + j, prec);
                acb_mul_mid(y, y, t, prec);
            }
        }

        mag_zero(arb_radref(acb_realref(y)));
        mag_zero(arb_radref(acb_imagref(y)));

        /* In the extremely unlikely event that y = 0, don't divide by 0. */
        if (arf_is_zero(arb_midref(acb_realref(y))) && arf_is_zero(arb_midref(acb_imagref(y))))
        {
            arf_set_ui_2exp_si(arb_midref(acb_realref(y)), 1, -prec);
            arf_set_ui_2exp_si(arb_midref(acb_imagref(y)), 1, -prec);
        }

        acb_inv_mid(t, y, prec);
        acb_mul_mid(t, t, x, prec);
        acb_sub_mid(roots + i, roots + i, t, prec);

        arf_get_mag(arb_radref(acb_realref(roots + i)), arb_midref(acb_realref(t)));
        arf_get_mag(arb_radref(acb_imagref(roots + i)), arb_midref(acb_imagref(t)));
    }

    acb_clear(x);
    acb_clear(y);
    acb_clear(t);
}

