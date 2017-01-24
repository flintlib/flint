/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

slong _fmpz_mpoly_sub1(fmpz * poly1, ulong * exps1,
                 const fmpz * poly2, const ulong * exps2, slong len2,
                           const fmpz * poly3, const ulong * exps3, slong len3)
{
   slong i = 0, j = 0, k = 0;

   while (i < len2 && j < len3)
   {
      if (exps2[i] < exps3[j])
      {
         fmpz_set(poly1 + k, poly2 + i);
         exps1[k] = exps2[i];
         i++;
      } else if (exps2[i] == exps3[j])
      {
         fmpz_sub(poly1 + k, poly2 + i, poly3 + j);
         exps1[k] = exps2[i];
         if (fmpz_is_zero(poly1 + k))
            k--;
         i++;
         j++;
      } else
      {
         fmpz_neg(poly1 + k, poly3 + j);
         exps1[k] = exps3[j];
         j++;         
      }
      k++;
   }

   while (i < len2)
   {
      fmpz_set(poly1 + k, poly2 + i);
      exps1[k] = exps2[i];
      i++;
      k++;
   }

   while (j < len3)
   {
      fmpz_neg(poly1 + k, poly3 + j);
      exps1[k] = exps3[j];
      j++;
      k++;
   }

   return k;
}

slong _fmpz_mpoly_sub(fmpz * poly1, ulong * exps1,
                  const fmpz * poly2, const ulong * exps2, slong len2,
                  const fmpz * poly3, const ulong * exps3, slong len3, slong N)
{
   slong i = 0, j = 0, k = 0;

   while (i < len2 && j < len3)
   {
      int cmp = mpoly_monomial_cmp(exps2 + i*N, exps3 + j*N, N);

      if (cmp < 0)
      {
         fmpz_set(poly1 + k, poly2 + i);
         mpoly_monomial_set(exps1 + k*N, exps2 + i*N, N);
         i++;
      } else if (cmp == 0)
      {
         fmpz_sub(poly1 + k, poly2 + i, poly3 + j);
         mpoly_monomial_set(exps1 + k*N, exps2 + i*N, N);
         if (fmpz_is_zero(poly1 + k))
            k--;
         i++;
         j++;
      } else
      {
         fmpz_neg(poly1 + k, poly3 + j);
         mpoly_monomial_set(exps1 + k*N, exps3 + j*N, N);
         j++;         
      }
      k++;
   }

   while (i < len2)
   {
      fmpz_set(poly1 + k, poly2 + i);
      mpoly_monomial_set(exps1 + k*N, exps2 + i*N, N);
      i++;
      k++;
   }

   while (j < len3)
   {
      fmpz_neg(poly1 + k, poly3 + j);
      mpoly_monomial_set(exps1 + k*N, exps3 + j*N, N);
      j++;
      k++;
   }

   return k;
}

void fmpz_mpoly_sub(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong len = 0, max_bits, N;
   ulong * ptr1, * ptr2;
   int free2, free3;

   max_bits = FLINT_MAX(poly2->bits, poly3->bits);
   N = (max_bits*ctx->n - 1)/FLINT_BITS + 1;
   
   ptr1 = _fmpz_mpoly_unpack_monomials(max_bits, poly2->exps, 
                                           poly2->bits, ctx->n, poly2->length);

   free2 = ptr1 != poly2->exps;

   ptr2 = _fmpz_mpoly_unpack_monomials(max_bits, poly3->exps, 
                                           poly3->bits, ctx->n, poly3->length);

   free3 = ptr2 != poly3->exps;

   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length, ctx);
      fmpz_mpoly_fit_bits(temp, max_bits, ctx);

      len = _fmpz_mpoly_sub(temp->coeffs, temp->exps, 
                    poly2->coeffs, ptr1, poly2->length,
                    poly3->coeffs, ptr2, poly3->length, N);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length, ctx);
      fmpz_mpoly_fit_bits(poly1, max_bits, ctx);

      len = _fmpz_mpoly_sub(poly1->coeffs, poly1->exps, 
                       poly2->coeffs, ptr1, poly2->length,
                       poly3->coeffs, ptr2, poly3->length, N);
   }
      
   if (free2)
      flint_free(ptr1);

   if (free3)
      flint_free(ptr2);

   _fmpz_mpoly_set_length(poly1, len, ctx);
}
