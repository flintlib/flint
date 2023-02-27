/*
    Copyright (C) 2015, 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"
#include "fq_mat.h"

slong fq_mat_reduce_row(fq_mat_t A, slong * P, slong * L, 
                                         slong m, const fq_ctx_t ctx)
{
   slong n = A->c, i, j, r, res = -WORD(1);
   fmpz_poly_t h;

   fmpz_poly_init(h);

   for (i = 0; i < n; i++)
   {
      if (i != 0)
         fq_reduce(fq_mat_entry(A, m, i), ctx);

      if (!fq_is_zero(fq_mat_entry(A, m, i), ctx))
      {
         r = P[i];
         if (r != -WORD(1))
         {
            for (j = i + 1; j < L[r]; j++)
            {
               fmpz_poly_mul(h, fq_mat_entry(A, r, j), fq_mat_entry(A, m, i));
               fmpz_poly_sub(fq_mat_entry(A, m, j), fq_mat_entry(A, m, j), h);
            }
 
            fq_zero(fq_mat_entry(A, m, i), ctx);
         } else
         {
            fq_inv(h, fq_mat_entry(A, m, i), ctx);
            fq_one(fq_mat_entry(A, m, i), ctx);
           
            for (j = i + 1; j < L[m]; j++)
            {
               fq_reduce(fq_mat_entry(A, m, j), ctx);

               fq_mul(fq_mat_entry(A, m, j), fq_mat_entry(A, m, j), h, ctx);
            }

            P[i] = m;

            res = i;

            break;
         }
      }
   }

   fmpz_poly_clear(h);
   
   return res;
}

