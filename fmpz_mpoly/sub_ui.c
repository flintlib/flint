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

void fmpz_mpoly_sub_ui(fmpz_mpoly_t poly1,
                 const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx)
{
   slong i, m = (poly2->bits*ctx->n - 1)/FLINT_BITS + 1;

   if (poly2->length == 0)
   {
      fmpz_mpoly_set_ui(poly1, c, ctx);

      if (c != 0)
         fmpz_neg(poly1->coeffs + 0, poly1->coeffs + 0);

      return;
   }

   if (c != 0)
   {
      if (mpoly_monomial_is_zero(poly2->exps + 0, m))
      {
         if (poly1 != poly2)
         {
            fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
            fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);

            for (i = 1; i < poly2->length; i++)
               fmpz_set(poly1->coeffs + i, poly2->coeffs + i);

            for (i = 0; i < poly2->length*m; i++)
               poly1->exps[i] = poly2->exps[i];

            _fmpz_mpoly_set_length(poly1, poly2->length, ctx);

            poly1->bits = poly2->bits;
         }

         fmpz_sub_ui(poly1->coeffs + 0, poly2->coeffs + 0, c);

         _fmpz_mpoly_renormalise(poly1, ctx);
      } else
      {
         fmpz_mpoly_fit_length(poly1, poly2->length + 1, ctx);

         if (poly1 != poly2)
         {
            fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);

            for (i = poly2->length - 1; i >= 0; i--)
               fmpz_set(poly1->coeffs + i + 1, poly2->coeffs + i);

            poly1->bits = poly2->bits;
         } else
         {
            for (i = poly2->length - 1; i >= 0; i--)
               poly1->coeffs[i + 1] = poly2->coeffs[i];
         }

         for (i = poly2->length*m - 1; i >= 0; i--)
            poly1->exps[i + m] = poly2->exps[i];
            
         _fmpz_mpoly_set_length(poly1, poly2->length + 1, ctx);

         for (i = 0; i < m; i++)
            poly1->exps[i] = 0;

         fmpz_set_ui(poly1->coeffs + 0, c);
         fmpz_neg(poly1->coeffs + 0, poly1->coeffs + 0);
      }
   }
}
