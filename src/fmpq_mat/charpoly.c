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

void _fmpq_mat_charpoly(fmpz * coeffs, fmpz_t den, const fmpq_mat_t mat)
{
   fmpz_t p;
   fmpz_mat_t zmat;
   slong n = mat->r, i;

   if (n == 0)
   {
      fmpz_one(coeffs + 0);
      fmpz_one(den);
      return;
   }

   fmpz_init(p);
   fmpz_mat_init(zmat, mat->r, mat->c);

   fmpq_mat_get_fmpz_mat_matwise(zmat, p, mat);

   _fmpz_mat_charpoly(coeffs, zmat);

   if (n > 1)
   {
      fmpz_mul(coeffs + 1, coeffs + 1, p);
      fmpz_mul(den, p, p);
   } else
      fmpz_set(den, p);

   for (i = 2; i < n; i++)
   {
      fmpz_mul(coeffs + i, coeffs + i, den);
      fmpz_mul(den, den, p);
   }

   fmpz_mul(coeffs + n, coeffs + n, den);

   fmpz_mat_clear(zmat);
   fmpz_clear(p);
}

void fmpq_mat_charpoly(fmpq_poly_t pol, const fmpq_mat_t mat)
{
   slong n = mat->r;

   if (mat->r != mat->c)
   {
       flint_throw(FLINT_ERROR, "Exception (fmpq_mat_charpoly).  Non-square matrix.\n");
   }

   fmpq_poly_fit_length(pol, n + 1);
   _fmpq_poly_set_length(pol, n + 1);

   _fmpq_mat_charpoly(pol->coeffs, pol->den, mat);

   fmpq_poly_canonicalise(pol);
}
