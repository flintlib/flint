/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A,
    slong len, const arb_t c, slong prec)
{
    arb_t r, t;
    slong i;

    if (len < 2)
    {
        arb_zero(R);
        return;
    }

    arb_init(r);
    arb_init(t);

    arb_set(t, A + len - 2);
    arb_set(Q + len - 2, A + len - 1);
    arb_set(r, Q + len - 2);

    /* TODO: avoid the extra assignments (but still support aliasing)  */
    for (i = len - 2; i > 0; i--)
    {
        arb_mul(r, r, c, prec);
        arb_add(r, r, t, prec);
        arb_set(t, A + i - 1);
        arb_set(Q + i - 1, r);
    }

    arb_mul(r, r, c, prec);
    arb_add(R, r, t, prec);

    arb_clear(r);
    arb_clear(t);
}

