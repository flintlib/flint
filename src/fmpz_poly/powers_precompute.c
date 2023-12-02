/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"

fmpz **
_fmpz_poly_powers_precompute(const fmpz * B, slong len)
{
   slong i;
   fmpz ** powers = flint_malloc(sizeof(fmpz *)*(2*len - 1));
   fmpz_poly_t pow, p;

   fmpz_poly_init2(pow, len);
   fmpz_poly_one(pow);
   fmpz_poly_init2(p, len - 1);

   for (i = 0; i < 2*len - 1; i++)
   {
      powers[i] = _fmpz_vec_init(len - 1);

      if (pow->length == len) /* reduce pow mod B */
      {
         _fmpz_vec_scalar_mul_fmpz(p->coeffs, B, len - 1, pow->coeffs + pow->length - 1);
         _fmpz_poly_set_length(p, len - 1);
         _fmpz_poly_normalise(p);
         fmpz_poly_sub(pow, pow, p);
         _fmpz_poly_set_length(pow, len - 1);
         _fmpz_poly_normalise(pow);
      }

      _fmpz_vec_set(powers[i], pow->coeffs, len - 1);
      fmpz_poly_shift_left(pow, pow, 1);
   }

   fmpz_poly_clear(pow);
   fmpz_poly_clear(p);

   return powers;
}

void fmpz_poly_powers_precompute(fmpz_poly_powers_precomp_t pinv,
                                                          fmpz_poly_t poly)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_powers_precompute). Division by zero.\n");
    }

    pinv->powers = _fmpz_poly_powers_precompute(poly->coeffs, poly->length);
    pinv->len = poly->length;
}
