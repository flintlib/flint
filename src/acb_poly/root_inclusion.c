/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

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

    _acb_poly_evaluate(t, poly, len, r, prec);
    acb_get_abs_ubound_arf(u, t, MAG_BITS);

    /* it could happen that we have an exact root, in which case
       we should avoid dividing by the derivative */
    if (!arf_is_zero(u))
    {
        _acb_poly_evaluate(t, polyder, len - 1, r, prec);
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

