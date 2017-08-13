/*
    Copyright (C) 2017 William Hart

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

void _fmpz_mpoly_reverse1(fmpz * poly1, ulong * exp1, const fmpz * poly2,
                                                 const ulong * exp2, slong len)
{
   slong i;
   ulong t;

   if (poly1 == poly2)
   {
      for (i = 0; i < len/2; i++)
      {
         fmpz_swap(poly1 + i, poly1 + len - i - 1);
         t = exp1[i];
         exp1[i] = exp1[len - i - 1];
         exp1[len - i - 1] = t;
      }
   } else
   {
      for (i = 0; i < len; i++)
      {
         fmpz_set(poly1 + i, poly2 + len - i - 1);
         exp1[i] = exp2[len - i - 1];
      }
   }
}
 
void _fmpz_mpoly_reverse(fmpz * poly1, ulong * exp1, const fmpz * poly2,
                                        const ulong * exp2, slong len, slong N)
{
   slong i;

   if (N == 1)
   {
      _fmpz_mpoly_reverse1(poly1, exp1, poly2, exp2, len);

      return;
   }

   if (poly1 == poly2)
   {
      for (i = 0; i < len/2; i++)
      {
         fmpz_swap(poly1 + i, poly1 + len - i - 1);
         mpoly_monomial_swap(exp1 + i*N, exp1 + (len - i - 1)*N, N);
      }
   } else
   {
      for (i = 0; i < len; i++)
      {
         fmpz_set(poly1 + i, poly2 + len - i - 1);
         mpoly_monomial_set(exp1 + i*N, exp2 + (len - i - 1)*N, N);
      }
   }
}
 
void fmpz_mpoly_reverse(fmpz_mpoly_t poly1,
                                fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
   slong N = (poly2->bits*ctx->n - 1)/FLINT_BITS + 1;

   if (poly1 != poly2)
   {
      fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
      fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);
   }

   _fmpz_mpoly_reverse(poly1->coeffs, poly1->exps, poly2->coeffs,
                                                poly2->exps, poly2->length, N);

   _fmpz_mpoly_set_length(poly1, poly2->length, ctx);
}
