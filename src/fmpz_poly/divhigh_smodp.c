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

void fmpz_poly_divhigh_smodp(fmpz * res, const fmpz_poly_t f,
                                  const fmpz_poly_t g, const fmpz_t p, slong n)
{
   fmpz_t d, cinv;
   slong i = 0, k, start = 0, len_g = g->length;
   fmpz_poly_t tf;

   fmpz_init(d);
   fmpz_init(cinv);

   fmpz_poly_init2(tf, f->length);

   fmpz_poly_set(tf, f);

   fmpz_gcdinv(d, cinv, g->coeffs + len_g - 1, p);
   if (!fmpz_is_one(d))
   {
      flint_throw(FLINT_ERROR, "Exception (fmpz_poly_divhigh_smodp). Impossible inverse.\n");
   }

   for (k = n - 1, i = f->length - len_g; k >= 0; i--, k--)
   {
      if (i < f->length - n)
         start++;

      fmpz_mul(res + k, tf->coeffs + i + len_g - 1, cinv);

      fmpz_smod(res + k, res + k, p);

      _fmpz_vec_scalar_submul_fmpz(tf->coeffs + i + start,
                                    g->coeffs + start, len_g - start, res + k);
      _fmpz_vec_scalar_smod_fmpz(tf->coeffs + i + start,
                                     tf->coeffs + i + start, len_g - start, p);
   }

   fmpz_poly_clear(tf);
   fmpz_clear(cinv);
   fmpz_clear(d);
}
