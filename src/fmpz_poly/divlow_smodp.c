/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void fmpz_poly_divlow_smodp(fmpz * res, const fmpz_poly_t f,
                                  const fmpz_poly_t g, const fmpz_t p, slong n)
{
   fmpz_t d, cinv;
   slong i = 0, k, zeroes;
   fmpz_poly_t tf;

   fmpz_init(d);
   fmpz_init(cinv);

   while (fmpz_is_zero(g->coeffs + i))
      i++;

   zeroes = i;

   fmpz_poly_init2(tf, n + zeroes);

   fmpz_poly_set(tf, f);

   if (fmpz_sgn(g->coeffs + zeroes) >= 0)
      fmpz_gcdinv(d, cinv, g->coeffs + zeroes, p);
   else
   {
      fmpz_t temp;

      fmpz_init(temp);

      fmpz_add(temp, g->coeffs + zeroes, p);
      fmpz_gcdinv(d, cinv, temp, p);

      fmpz_clear(temp);
   }

   if (!fmpz_is_one(d))
   {
      flint_throw(FLINT_ERROR, "Exception (fmpz_poly_divlow_smodp). Impossible inverse.\n");
   }

   for (k = 0; k < n; i++, k++)
   {
      fmpz_mul(res + k, tf->coeffs + i, cinv);

      fmpz_smod(res + k, res + k, p);

      _fmpz_vec_scalar_submul_fmpz(tf->coeffs + i,  g->coeffs + zeroes,
                                FLINT_MIN(g->length - zeroes, n - k), res + k);
      _fmpz_vec_scalar_smod_fmpz(tf->coeffs + i, tf->coeffs + i,
                                      FLINT_MIN(g->length - zeroes, n - k), p);
   }

   fmpz_poly_clear(tf);
   fmpz_clear(cinv);
   fmpz_clear(d);
}
