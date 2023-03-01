/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
acb_poly_clear(acb_poly_t poly)
{
    slong i;

    for (i = 0; i < poly->alloc; i++)
        acb_clear(poly->coeffs + i);

    flint_free(poly->coeffs);
}
