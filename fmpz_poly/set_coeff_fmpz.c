/*
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
#include "fmpz_poly.h"

void
fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, slong n, const fmpz_t x)
{
    if (fmpz_is_zero(x))
    {
       if (n >= poly->length)
          return;

       fmpz_zero(poly->coeffs + n);

       if (n == poly->length - 1) /* only necessary when setting leading coefficient */
          _fmpz_poly_normalise(poly);
    }
    else
    {
        fmpz_poly_fit_length(poly, n + 1);

        if (n + 1 > poly->length)
        {
           slong i;
           
           for (i = poly->length; i < n; i++)
               fmpz_zero(poly->coeffs + i);
           
           poly->length = n + 1;
        }

        fmpz_set(poly->coeffs + n, x);
    }
}
