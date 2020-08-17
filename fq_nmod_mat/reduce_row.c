/*
    Copyright (C) 2015, 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fq_nmod.h"
#include "fq_nmod_mat.h"

slong fq_nmod_mat_reduce_row_KS(fq_nmod_mat_t A, slong * P, slong * L,
                                         slong m, const fq_nmod_ctx_t ctx)
{
   slong n = A->c, i, j, r, bits, res = -WORD(1);
   fq_nmod_t h;
   fmpz * mvec;
   fmpz_t mz, rz;

   bits = FLINT_BIT_COUNT(ctx->p)*2 + FLINT_BIT_COUNT(m + 1) +
          FLINT_BIT_COUNT(fq_nmod_ctx_degree(ctx) + 1);

   fq_nmod_init(h, ctx);
   fmpz_init(mz);
   fmpz_init(rz);
   mvec = (fmpz *) _fmpz_vec_init(n);

   for (i = 0; i < n; i++)
      fq_nmod_bit_pack(mvec + i, fq_nmod_mat_entry(A, m, i), bits, ctx);

   for (i = 0; i < n; i++)
   {
      if (i != 0)
         fq_nmod_bit_unpack(fq_nmod_mat_entry(A, m, i),
                                            mvec + i, bits, ctx);

      if (!fq_nmod_is_zero(fq_nmod_mat_entry(A, m, i), ctx))
      {
         r = P[i];
         if (r != -WORD(1))
         {
            fq_nmod_neg(h, fq_nmod_mat_entry(A, m, i), ctx);
            fq_nmod_bit_pack(mz, h, bits, ctx);

            for (j = i + 1; j < L[r]; j++)
            {
               fq_nmod_bit_pack(rz, fq_nmod_mat_entry(A, r, j), bits, ctx);
               fmpz_mul(rz, rz, mz);

               fmpz_add(mvec + j, mvec + j, rz);
            }

            fq_nmod_zero(fq_nmod_mat_entry(A, m, i), ctx);
         } else
         {
            fq_nmod_inv(h, fq_nmod_mat_entry(A, m, i), ctx);
            fq_nmod_one(fq_nmod_mat_entry(A, m, i), ctx);

            for (j = i + 1; j < L[m]; j++)
            {
               fq_nmod_bit_unpack(fq_nmod_mat_entry(A, m, j),
                                            mvec + j, bits, ctx);

               fq_nmod_mul(fq_nmod_mat_entry(A, m, j), fq_nmod_mat_entry(A, m, j), h, ctx
);
            }

            P[i] = m;

            res = i;

            break;
         }
      }
   }

   fq_nmod_clear(h, ctx);
   fmpz_clear(mz);
   fmpz_clear(rz);
   _fmpz_vec_clear(mvec, n);

   return res;
}

slong fq_nmod_mat_reduce_row(fq_nmod_mat_t A, slong * P, slong * L, 
                                         slong m, const fq_nmod_ctx_t ctx)
{
   slong n = A->c, i, j, r;
   nmod_poly_t h;

   if (m > 10 && fq_nmod_ctx_degree(ctx) > 6)
      return fq_nmod_mat_reduce_row_KS(A, P, L, m, ctx);

   nmod_poly_init(h, ctx->p);

   for (i = 0; i < n; i++)
   {
      if (i != 0)
         fq_nmod_reduce(fq_nmod_mat_entry(A, m, i), ctx);

      if (!fq_nmod_is_zero(fq_nmod_mat_entry(A, m, i), ctx))
      {
         r = P[i];
         if (r != -WORD(1))
         {
            for (j = i + 1; j < L[r]; j++)
            {
               nmod_poly_mul(h, fq_nmod_mat_entry(A, r, j), fq_nmod_mat_entry(A, m, i));
               nmod_poly_sub(fq_nmod_mat_entry(A, m, j), fq_nmod_mat_entry(A, m, j), h);
            }
 
            fq_nmod_zero(fq_nmod_mat_entry(A, m, i), ctx);
         } else
         {
            fq_nmod_inv(h, fq_nmod_mat_entry(A, m, i), ctx);
            fq_nmod_one(fq_nmod_mat_entry(A, m, i), ctx);
           
            for (j = i + 1; j < L[m]; j++)
            {
               fq_nmod_reduce(fq_nmod_mat_entry(A, m, j), ctx);

               fq_nmod_mul(fq_nmod_mat_entry(A, m, j), fq_nmod_mat_entry(A, m, j), h, ctx);
            }

            P[i] = m;

            nmod_poly_clear(h);
       
            return i;
         }
      }
   }

   for (j = i + 1; j < L[m]; j++)
      fq_nmod_reduce(fq_nmod_mat_entry(A, m, j), ctx);

   nmod_poly_clear(h);
   
   return -WORD(1);
}

