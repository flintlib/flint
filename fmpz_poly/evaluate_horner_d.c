/*
    Copyright (C) 2016 William Hart

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

double _fmpz_poly_evaluate_horner_d(const fmpz * poly, slong n, double d)
{
   slong i;
   double ans;

   ans = fmpz_get_d(poly + n - 1);
   
   for (i = n - 2; i >= 0L; i--)
   {
      ans *= d;
      ans += fmpz_get_d(poly + i);
   }

   return ans;
}

double fmpz_poly_evaluate_horner_d(const fmpz_poly_t poly, double d)
{
   if (poly->length == 0)
      return 0.0;

   return _fmpz_poly_evaluate_horner_d(poly->coeffs, poly->length, d);
}
