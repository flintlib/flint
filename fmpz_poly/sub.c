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
_fmpz_poly_sub(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2,
               slong len2)
{
    slong i, min = FLINT_MIN(len1, len2);

    for (i = 0; i < min; i++)   /* subtract up to the length of the shorter poly */
        fmpz_sub(res + i, poly1 + i, poly2 + i);

    if (poly1 != res)           /* copy any remaining coefficients from poly1 */
        for (i = min; i < len1; i++)
            fmpz_set(res + i, poly1 + i);

    /* careful, it is *always* necessary to negate coeffs from poly2, even if this is already res */
    for (i = min; i < len2; i++)
        fmpz_neg(res + i, poly2 + i);
}

void
fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    fmpz_poly_fit_length(res, max);

    _fmpz_poly_sub(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length);

    _fmpz_poly_set_length(res, max);
    _fmpz_poly_normalise(res);  /* there may have been cancellation */
}

void fmpz_poly_sub_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c)
{
   if (poly->length == 0)
   {
      fmpz_poly_set_si(res, c);
      fmpz_poly_neg(res, res);
   }
   else
   {
      fmpz_poly_set(res, poly);

      if (c < 0)
         fmpz_add_ui(res->coeffs + 0, res->coeffs + 0, -c);
      else
         fmpz_sub_ui(res->coeffs + 0, res->coeffs + 0, c);

      _fmpz_poly_normalise(res);
   }
}

void fmpz_poly_si_sub(fmpz_poly_t res, slong c, const fmpz_poly_t poly)
{
    if (poly->length == 0)
        fmpz_poly_set_si(res, c);
    else
    {
        fmpz_poly_neg(res, poly);

        if (c < 0)
            fmpz_sub_ui(res->coeffs + 0, res->coeffs + 0, -c);
         else
            fmpz_add_ui(res->coeffs + 0, res->coeffs + 0, c);

        _fmpz_poly_normalise(res);
    }
}

void fmpz_poly_sub_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c)
{
    if (poly->length == 0)
    {
        fmpz_poly_set_fmpz(res, c);
        fmpz_neg(res->coeffs + 0, res->coeffs + 0);
    }
    else
    {
        fmpz_poly_set(res, poly);
        fmpz_sub(res->coeffs + 0, res->coeffs + 0, c);
        _fmpz_poly_normalise(res);
    }
}

void fmpz_poly_fmpz_sub(fmpz_poly_t res, fmpz_t c, const fmpz_poly_t poly)
{
    if (poly->length == 0)
        fmpz_poly_set_fmpz(res, c);
    else
    {
        fmpz_poly_neg(res, poly);
        fmpz_add(res->coeffs + 0, res->coeffs + 0, c);
        _fmpz_poly_normalise(res);
    }
}
