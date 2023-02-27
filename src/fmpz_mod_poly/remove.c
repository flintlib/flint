/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

ulong fmpz_mod_poly_remove(fmpz_mod_poly_t f, const fmpz_mod_poly_t g,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_t q, r;
    ulong i = 0;

    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_init(r, ctx);

    while (1)
    {
        if (f->length < g->length)
            break;
        fmpz_mod_poly_divrem(q, r, f, g, ctx);
        if (r->length == 0)
            fmpz_mod_poly_swap(q, f, ctx);
        else
            break;
        i++;
    }

    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_poly_clear(r, ctx);

    return i;
}
