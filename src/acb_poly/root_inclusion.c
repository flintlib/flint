/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

/* Horner with a tighter error bound based on a disk representation instead of
 * the pair of a rectangle. This is critical for evaluating high degree
 * polynomials on complex points */
static void _acb_poly_evaluate_stable(acb_t y, acb_srcptr poly, slong len, const acb_t x, slong prec)
{
    mag_t ax, e, f;
    slong i;

    mag_init(ax);
    mag_init(e);
    mag_init(f);

    acb_get_mag(ax, x);
    acb_zero(y);

    if(len>0) {
        acb_set(y, poly + len - 1);
        mag_hypot(e, arb_radref(acb_realref(y)),
                     arb_radref(acb_imagref(y)));

        for( i = len - 2; i >= 0; i--) {
            acb_mul(y, y, x, prec);
            acb_add(y, y, poly + i, prec);
            mag_mul(e, e, ax);
            /* crude upper bound on hypot here because it doesn't accumulate */
            mag_add(f, arb_radref(acb_realref(y)),
                       arb_radref(acb_imagref(y)));
            mag_add(e, e, f);
            acb_get_mid(y, y);
        }

        mag_set(arb_radref(acb_realref(y)), e);
        mag_set(arb_radref(acb_imagref(y)), e);
    }

    mag_clear(f);
    mag_clear(e);
    mag_clear(ax);
}


void
_acb_poly_root_inclusion(acb_t r, const acb_t m,
    acb_srcptr poly,
    acb_srcptr polyder, slong len, slong prec)
{
    acb_t t;
    arf_t u, v;

    acb_init(t);
    arf_init(u);
    arf_init(v);

    acb_set(r, m);
    mag_zero(arb_radref(acb_realref(r)));
    mag_zero(arb_radref(acb_imagref(r)));

    _acb_poly_evaluate_stable(t, poly, len, r, prec);
    acb_get_abs_ubound_arf(u, t, MAG_BITS);

    /* it could happen that we have an exact root, in which case
       we should avoid dividing by the derivative */
    if (!arf_is_zero(u))
    {
        _acb_poly_evaluate_stable(t, polyder, len - 1, r, prec);
        acb_inv(t, t, MAG_BITS);
        acb_get_abs_ubound_arf(v, t, MAG_BITS);

        arf_mul(u, u, v, MAG_BITS, ARF_RND_UP);
        arf_mul_ui(u, u, len - 1, MAG_BITS, ARF_RND_UP);
    }

    arf_get_mag(arb_radref(acb_realref(r)), u);
    arf_get_mag(arb_radref(acb_imagref(r)), u);

    arf_clear(u);
    arf_clear(v);
    acb_clear(t);
}
