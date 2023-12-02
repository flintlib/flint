/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"
#include "fmpq_mat.h"

slong _fmpq_mat_minpoly(fmpz * coeffs, fmpz_t den, const fmpq_mat_t mat)
{
   fmpz_t p;
   fmpz_mat_t zmat;
   slong n = mat->r, len, i;

   if (n == 0)
   {
      fmpz_one(coeffs + 0);
      fmpz_one(den);
      return 1;
   }

   fmpz_init(p);
   fmpz_mat_init(zmat, mat->r, mat->c);

   fmpq_mat_get_fmpz_mat_matwise(zmat, p, mat);

   len = _fmpz_mat_minpoly(coeffs, zmat);

   if (len > 2)
   {
      fmpz_mul(coeffs + 1, coeffs + 1, p);
      fmpz_mul(den, p, p);
   } else
      fmpz_set(den, p);

   for (i = 2; i < len - 1; i++)
   {
      fmpz_mul(coeffs + i, coeffs + i, den);
      fmpz_mul(den, den, p);
   }

   fmpz_mul(coeffs + len - 1, coeffs + len - 1, den);

   fmpz_mat_clear(zmat);
   fmpz_clear(p);

   return len;
}

void fmpq_mat_minpoly(fmpq_poly_t pol, const fmpq_mat_t mat)
{
   slong n = mat->r, len;

   if (mat->r != mat->c)
   {
       flint_throw(FLINT_ERROR, "Exception (fmpq_mat_minpoly).  Non-square matrix.\n");
   }

   fmpq_poly_fit_length(pol, n + 1);

   len = _fmpq_mat_minpoly(pol->coeffs, pol->den, mat);

   _fmpq_poly_set_length(pol, len);
   fmpq_poly_canonicalise(pol);
}
