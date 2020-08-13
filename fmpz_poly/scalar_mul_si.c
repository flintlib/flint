/*
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
#include "fmpz_poly.h"

void
fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, const fmpz_poly_t poly2, slong x)
{
    slong i;

    /* Either scalar or input poly is zero */
    if ((x == WORD(0)) || (poly2->length == 0))
    {
        fmpz_poly_zero(poly1);
        return;
    }

    /* Special case, multiply by 1 */
    if (x == WORD(1))
    {
        fmpz_poly_set(poly1, poly2);
        return;
    }

    /* Special case, multiply by -1 */
    if (x == WORD(-1))
    {
        fmpz_poly_neg(poly1, poly2);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);

    for (i = 0; i < poly2->length; i++)
        fmpz_mul_si(poly1->coeffs + i, poly2->coeffs + i, x);

    _fmpz_poly_set_length(poly1, poly2->length);
}
