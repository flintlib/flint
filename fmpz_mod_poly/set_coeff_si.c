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
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_set_coeff_si(fmpz_mod_poly_t poly, slong n, slong x,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        flint_mpn_zero((mp_ptr) (poly->coeffs + poly->length), n - poly->length);
        poly->length = n + 1;
    }

    fmpz_set_si(poly->coeffs + n, x);
    fmpz_mod(poly->coeffs + n, poly->coeffs + n, fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_normalise(poly);
}

