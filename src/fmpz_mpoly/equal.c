/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

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
#include "fmpz_mpoly.h"

int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1,
                     const fmpz * poly2, const ulong * exps2, slong n, slong N)
{
   slong i;

   if (poly1 != poly2)
   {
      for (i = 0; i < n; i++)
      {
         if (!fmpz_equal(poly1 + i, poly2 + i))
            return 0;
      }
   }

   if (exps1 != exps2)
   {
      for (i = 0; i < n*N; i++)
      {
         if (exps1[i] != exps2[i])
            return 0;
      }
   }

   return 1;
}

int fmpz_mpoly_equal(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   ulong * ptr1 = poly1->exps, * ptr2 = poly2->exps;
   slong max_bits, N;
   int r, free1 = 0, free2 = 0;

   if (poly1 == poly2)
      return 1;

   if (poly1->length != poly2->length)
      return 0;

   max_bits = FLINT_MAX(poly1->bits, poly2->bits);
   N = mpoly_words_per_exp(max_bits, ctx->minfo);

   if (max_bits > poly1->bits)
   {
      free1 = 1;
      ptr1 = (ulong *) flint_malloc(N*poly1->length*sizeof(ulong));
      mpoly_repack_monomials(ptr1, max_bits, poly1->exps, poly1->bits,
                                                    poly1->length, ctx->minfo);
   }

   if (max_bits > poly2->bits)
   {
      free2 = 1;
      ptr2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_repack_monomials(ptr2, max_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
   }

   r = _fmpz_mpoly_equal(poly1->coeffs, ptr1,
                                        poly2->coeffs, ptr2, poly2->length, N);

   if (free1)
      flint_free(ptr1);

   if (free2)
      flint_free(ptr2);

   return r;
}
