/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"

static void
acb_poly_taylor_shift_si_trunc(acb_poly_t res, const acb_poly_t pol,
                               const slong n, slong ord, slong prec)
{
    /* todo: faster code, especially for small ord */
    acb_t a;
    acb_init(a);

    acb_set_si(a, n);
    acb_poly_taylor_shift(res, pol, a, prec);
    acb_poly_truncate(res, ord);

    acb_clear(a);
}

void
acb_ode_bound_rat_ref_vec(mag_ptr res,
                  const acb_poly_struct * num, slong len, const acb_poly_t ind,
                  slong n, slong mult, slong ord, slong prec)
{
    acb_poly_t inv, ser;
    mag_t m;

    if (n == 0)
    {
        for (slong i = 0; i < len; i++)
            mag_zero(res + i);
        return;
    }

    acb_poly_init(inv);
    acb_poly_init(ser);
    mag_init(m);

    acb_poly_taylor_shift_si_trunc(inv, ind, n, ord + mult, prec);
    acb_poly_shift_right(inv, inv, mult);
    acb_poly_inv_series(inv, inv, ord, prec);

    for (slong i = 0; i < len; i++)
    {
        acb_poly_taylor_shift_si_trunc(ser, num + i, n, ord, prec);
        acb_poly_mullow(ser, ser, inv, ord, prec);

        mag_zero(res + i);
        for (slong j = 0; j < ser->length; j++)
        {
            acb_get_mag(m, ser->coeffs + j);
            mag_add(res + i, res + i, m);
        }
        mag_mul_ui(res + i, res + i, n);
    }

    acb_poly_clear(inv);
    acb_poly_clear(ser);
    mag_clear(m);
}
