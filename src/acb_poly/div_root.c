/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_div_root(acb_ptr Q, acb_t R, acb_srcptr A,
    slong len, const acb_t c, slong prec)
{
    acb_t r, t;
    slong i;

    if (len < 2)
    {
        acb_zero(R);
        return;
    }

    acb_init(r);
    acb_init(t);

    acb_set(t, A + len - 2);
    acb_set(Q + len - 2, A + len - 1);
    acb_set(r, Q + len - 2);

    /* TODO: avoid the extra assignments (but still support aliasing)  */
    for (i = len - 2; i > 0; i--)
    {
        acb_mul(r, r, c, prec);
        acb_add(r, r, t, prec);
        acb_set(t, A + i - 1);
        acb_set(Q + i - 1, r);
    }

    acb_mul(r, r, c, prec);
    acb_add(R, r, t, prec);

    acb_clear(r);
    acb_clear(t);
}
