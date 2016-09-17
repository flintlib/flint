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

int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1,
                     const fmpz * poly2, const ulong * exps2, slong n, slong m)
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
      for (i = 0; i < n*m; i++)
      {
         if (exps1[i] != exps2[i])
            return 0;
      }
   }

   return 1;
}

int fmpz_mpoly_equal(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   ulong * ptr1, * ptr2;
   slong max_bits, m;
   int r;

   if (poly1 == poly2)
      return 1;

   if (poly1->length != poly2->length)
      return 0;

   max_bits = FLINT_MAX(poly1->bits, poly2->bits);
   m = (max_bits*ctx->n - 1)/FLINT_BITS + 1;

   ptr1 = _fmpz_mpoly_unpack_monomials(max_bits, poly1->exps, 
                                           poly1->bits, ctx->n, poly1->length);

   ptr2 = _fmpz_mpoly_unpack_monomials(max_bits, poly2->exps, 
                                           poly2->bits, ctx->n, poly2->length);

   r = _fmpz_mpoly_equal(poly1->coeffs, ptr1,
                                        poly2->coeffs, ptr2, poly2->length, m);

   if (ptr1 != poly1->exps)
      flint_free(ptr1);

   if (ptr2 != poly2->exps)
      flint_free(ptr2);

   return r;
}
