/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_set(fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (poly1 != poly2)         /* Aliasing is trivial */
    {
        slong i, len = poly2->length;

        fmpz_mod_poly_fit_length(poly1, len, ctx);

        for (i = 0; i < len; i++)
            fmpz_set(poly1->coeffs + i, poly2->coeffs + i);

        _fmpz_mod_poly_set_length(poly1, len);
    }
}
