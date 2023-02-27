/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
gr_poly_truncate(gr_poly_t poly, slong newlen, gr_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        _gr_poly_set_length(poly, newlen, ctx);
        _gr_poly_normalise(poly, ctx);
    }

    return GR_SUCCESS;
}
