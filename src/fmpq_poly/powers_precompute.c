/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

fmpq_poly_struct *
_fmpq_poly_powers_precompute(const fmpz * B, const fmpz_t denB, slong len)
{
   slong i;
   fmpq_poly_struct * powers = flint_malloc(sizeof(fmpq_poly_struct)*(2*len - 1));
   fmpq_poly_t pow, p;

   fmpq_poly_init2(pow, len);
   fmpq_poly_one(pow);
   fmpq_poly_init2(p, len - 1);

   for (i = 0; i < 2*len - 1; i++)
   {
      fmpq_poly_init(powers + i);

      if (pow->length == len) /* reduce pow mod B */
      {
         fmpz_mul(fmpq_poly_denref(p), B + len - 1, fmpq_poly_denref(pow));
         _fmpz_vec_scalar_mul_fmpz(fmpq_poly_numref(p), B, len - 1,
             fmpq_poly_numref(pow) + len - 1);
         _fmpq_poly_set_length(p, len - 1);
         _fmpq_poly_normalise(p);
         fmpq_poly_canonicalise(p);
         fmpq_poly_sub(pow, pow, p);
         _fmpq_poly_set_length(pow, len - 1);
         _fmpq_poly_normalise(pow);
         fmpq_poly_canonicalise(pow);
      }

      fmpq_poly_set(powers + i, pow);
      fmpq_poly_shift_left(pow, pow, 1);
   }

   fmpq_poly_clear(pow);
   fmpq_poly_clear(p);

   return powers;
}

void fmpq_poly_powers_precompute(fmpq_poly_powers_precomp_t pinv,
                                                          fmpq_poly_t poly)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_powers_precompute). Division by zero.\n");
    }

    pinv->powers = _fmpq_poly_powers_precompute(fmpq_poly_numref(poly),
                                         fmpq_poly_denref(poly), poly->length);
    pinv->len = poly->length;
}
