/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A)
{
    slong i, j, r, c;

    mag_t s, t;

    r = acb_mat_nrows(A);
    c = acb_mat_ncols(A);

    mag_zero(b);

    if (r == 0 || c == 0)
        return;

    mag_init(s);
    mag_init(t);

    for (i = 0; i < r; i++)
    {
        mag_zero(s);

        for (j = 0; j < c; j++)
        {
            acb_get_mag(t, acb_mat_entry(A, i, j));
            mag_add(s, s, t);
        }

        mag_max(b, b, s);
    }

    mag_clear(s);
    mag_clear(t);
}

