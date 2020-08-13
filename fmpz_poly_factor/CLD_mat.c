/*
    Copyright (C) 2011, 2016 William Hart
    Copyright (C) 2011 Andy Novocin

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
#include "fmpz_mat.h"

#include "fmpz_mod_poly.h"

slong _fmpz_poly_factor_CLD_mat(fmpz_mat_t res, const fmpz_poly_t f,
                              fmpz_poly_factor_t lifted_fac, fmpz_t P, ulong k)
{
   /*
      The product of the lifted factors g is equal to f modulo P. Here P = p^a.
      We will compute the top k and bottom k coefficients of each logarithmic
      derivative fg'/g. The results are stored in res, along with an extra
      row for the CLD bounds for that column. The matrix res is required to be
      initialised to be of size (r + 1, 2k).
   */

   slong i, zeroes, bound, lo_n, hi_n, r = lifted_fac->num;
   slong bit_r = FLINT_MAX(r, 20);
   fmpz_poly_t gd, gcld, temp;
   fmpz_poly_t trunc_f, trunc_fac; /* don't initialise trunc_f, trunc_fac */
   fmpz_t t;

   /* insert CLD bounds in last row of matrix */

   for (i = 0; i < k; i++)
   {
      fmpz_poly_CLD_bound(res->rows[r] + i, f, i);
      fmpz_poly_CLD_bound(res->rows[r] + 2*k - i - 1, f, f->length - i - 2);
   }

   /* we exclude columns in the middle for which CLD bounds are too large */

   fmpz_init(t);

   bound = fmpz_bits(P) - bit_r - bit_r/2; /* log_2(p^a / 2^{1.5r}) */

   for (lo_n = 0; lo_n < k; lo_n++)
   {
      fmpz_mul_ui(t, res->rows[r] + lo_n, (slong) sqrt(f->length));

      if (fmpz_bits(t) > bound)
         break;
   }

   fmpz_zero(t);

   for (hi_n = 0; hi_n < k; hi_n++)
   {
      fmpz_mul_ui(t, res->rows[r] + 2*k - hi_n - 1, (slong) sqrt(f->length));

      if (fmpz_bits(t) > bound)
         break;
   }

   fmpz_clear(t);

   /* now insert data into matrix */

   fmpz_poly_init(gd);
   fmpz_poly_init(gcld);
   /* do not initialise trunc_f */
   /* do not initialise trunc_fac */

   if (lo_n > 0)
   {
      for (i = 0; i < r; i++)
      {
         zeroes = 0;
         while (fmpz_is_zero(lifted_fac->p[i].coeffs + zeroes))
            zeroes++;

         fmpz_poly_attach_truncate(trunc_fac, lifted_fac->p + i, lo_n + zeroes + 1);
         fmpz_poly_derivative(gd, trunc_fac);
         fmpz_poly_mullow(gcld, f, gd, lo_n + zeroes);
         fmpz_poly_divlow_smodp(res->rows[i], gcld, trunc_fac, P, lo_n);
      }      
   }

   if (hi_n > 0)
   {
      fmpz_poly_init(temp);

      fmpz_poly_attach_shift(trunc_f, f, f->length - hi_n);
      
      for (i = 0; i < r; i++)
      {
         slong len = lifted_fac->p[i].length - hi_n - 1;

         if (len < 0)
         {
            fmpz_poly_shift_left(temp, lifted_fac->p + i, -len);
            fmpz_poly_derivative(gd, temp);
            fmpz_poly_mulhigh_n(gcld, trunc_f, gd, hi_n);
            fmpz_poly_divhigh_smodp(res->rows[i] + lo_n, gcld, temp, P, hi_n);
         } else
         {
            fmpz_poly_attach_shift(trunc_fac, lifted_fac->p + i, len);
            fmpz_poly_derivative(gd, trunc_fac);
            fmpz_poly_mulhigh_n(gcld, trunc_f, gd, hi_n);
            fmpz_poly_divhigh_smodp(res->rows[i] + lo_n, gcld, trunc_fac, P, hi_n);
         }
      }

      fmpz_poly_clear(temp);      
   }

   if (hi_n > 0)
   {
      /* move bounds into correct columns */
      for (i = 0; i < hi_n; i++)
         fmpz_set(res->rows[r] + lo_n + i, res->rows[r] + 2*k - hi_n + i);
   }

   /* do not clear trunc_fac */
   /* do not clear trunc_f */
   fmpz_poly_clear(gd);
   fmpz_poly_clear(gcld);

   return lo_n + hi_n;
}

