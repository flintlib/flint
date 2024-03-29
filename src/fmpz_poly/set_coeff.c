/*
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

void
fmpz_poly_set_coeff_si(fmpz_poly_t poly, slong n, slong x)
{
    if (x == 0)
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

        fmpz_set_si(poly->coeffs + n, x);
    }
}

void
fmpz_poly_set_coeff_ui(fmpz_poly_t poly, slong n, ulong x)
{
    if (x == 0)
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

        fmpz_set_ui(poly->coeffs + n, x);
    }
}
