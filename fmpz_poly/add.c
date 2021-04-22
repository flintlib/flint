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
_fmpz_poly_add(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2,
               slong len2)
{
    slong i, min = FLINT_MIN(len1, len2);

    for (i = 0; i < min; i++)   /* add up to the length of the shorter poly */
        fmpz_add(res + i, poly1 + i, poly2 + i);

    if (poly1 != res)           /* copy any remaining coefficients from poly1 */
        for (i = min; i < len1; i++)
            fmpz_set(res + i, poly1 + i);

    if (poly2 != res)           /* copy any remaining coefficients from poly2 */
        for (i = min; i < len2; i++)
            fmpz_set(res + i, poly2 + i);
}

void
fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    fmpz_poly_fit_length(res, max);

    _fmpz_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length);

    _fmpz_poly_set_length(res, max);
    _fmpz_poly_normalise(res);  /* there may have been cancellation */
}

void fmpz_poly_add_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c)
{
    if (poly->length == 0)
        fmpz_poly_set_si(res, c);
    else
    {
        fmpz_poly_set(res, poly);

        if (c < 0)
            fmpz_sub_ui(res->coeffs + 0, res->coeffs + 0, -c);
        else
            fmpz_add_ui(res->coeffs + 0, res->coeffs + 0, c);

        _fmpz_poly_normalise(res);
    }
}

void fmpz_poly_add_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c)
{
    if (poly->length == 0)
        fmpz_poly_set_fmpz(res, c);
    else
    {
        fmpz_poly_set(res, poly);
        fmpz_add(res->coeffs + 0, res->coeffs + 0, c);
        _fmpz_poly_normalise(res);
    }
}
