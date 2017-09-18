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

void fmpz_mpoly_add_fmpz(fmpz_mpoly_t poly1,
          const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
   slong i, N;
   slong len2 = poly2->length;

   if (len2 == 0)
   {
      fmpz_mpoly_set_fmpz(poly1, c, ctx);

      return;
   }

   if (!fmpz_is_zero(c))
   {
      N = words_per_exp(ctx->n, poly2->bits);

      if (mpoly_monomial_is_zero(poly2->exps + (len2 - 1)*N, N))
      {
         if (poly1 != poly2)
         {
            fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
            fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);

            for (i = 0; i < len2 - 1; i++)
               fmpz_set(poly1->coeffs + i, poly2->coeffs + i);

            for (i = 0; i < len2*N; i++)
               poly1->exps[i] = poly2->exps[i];

            _fmpz_mpoly_set_length(poly1, poly2->length, ctx);

            poly1->bits = poly2->bits;
         }

         fmpz_add(poly1->coeffs + len2 - 1, poly2->coeffs + len2 - 1, c);

         if (fmpz_is_zero(poly1->coeffs + len2 - 1))
            _fmpz_mpoly_set_length(poly1, len2 - 1, ctx);
      } else
      {
         fmpz_mpoly_fit_length(poly1, len2 + 1, ctx);

         if (poly1 != poly2)
         {
            fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);
            poly1->bits = poly2->bits;

            for (i = 0; i < len2; i++)
               fmpz_set(poly1->coeffs + i, poly2->coeffs + i);

            for (i = 0; i < len2*N; i++)
               poly1->exps[i] = poly2->exps[i];
         } 
            
         for (i = 0; i < N; i++)
            poly1->exps[len2*N + i] = 0;

         fmpz_set(poly1->coeffs + len2, c);

         _fmpz_mpoly_set_length(poly1, poly2->length + 1, ctx);
      }
   } else if (poly1 != poly2)
      fmpz_mpoly_set(poly1, poly2, ctx);
}
