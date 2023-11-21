/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec)
{
    acb_t s, t;
    mag_t err, C;

    acb_init(s);
    acb_init(t);
    mag_init(err);
    mag_init(C);

    if (n < 0)
        n = acb_hypgeom_pfq_choose_n(a, p, b, q, z, prec);

    acb_hypgeom_pfq_sum(s, t, a, p, b, q, z, n, prec);

    if (!acb_is_zero(t))
    {
        acb_hypgeom_pfq_bound_factor(C, a, p, b, q, z, n);
        acb_get_mag(err, t);
        mag_mul(err, err, C);

        if (_acb_vec_is_real(a, p) && _acb_vec_is_real(b, q) && acb_is_real(z))
            arb_add_error_mag(acb_realref(s), err);
        else
            acb_add_error_mag(s, err);
    }

    acb_set(res, s);

    acb_clear(s);
    acb_clear(t);
    mag_clear(err);
    mag_clear(C);
}

