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

slong fmpz_mpoly_get_term_si(const fmpz_mpoly_t poly,
                                 ulong const * exp, const fmpz_mpoly_ctx_t ctx)
{
   slong i, N, index, bits, exp_bits;
   int exists;
   ulong maskhi, masklo;
   ulong sum = 0, max_exp = 0;
   slong c;
   ulong * packed_exp;
   int deg, rev;

   TMP_INIT;

   TMP_START;

   degrev_from_ord(deg, rev, ctx->ord);

   if (deg)
   {
      for (i = 0; i < ctx->n - 1; i++)
      {  
         sum += exp[i];

         if (sum < exp[i])
            flint_throw(FLINT_EXPOF,
                              "Exponent overflow in fmpz_mpoly_get_term_fmpz");
      }

      max_exp = sum;
   } else
   {
      for (i = 0; i < ctx->n - 1; i++)
      {
         if (exp[i] > max_exp)
            max_exp = exp[i];
      }
   }
      
   if (0 > (slong) max_exp)
      flint_throw(FLINT_EXPOF,
                              "Exponent overflow in fmpz_mpoly_get_term_fmpz");

   /* compute number of bits to store maximum degree */
   bits = FLINT_BIT_COUNT(max_exp);

   exp_bits = 8;
   while (bits >= exp_bits) /* extra bit required for signs */
      exp_bits *= 2;

   if (exp_bits > poly->bits) /* exponent too large to be poly exponent */
      return 0;
   
   masks_from_bits_ord(maskhi, masklo, poly->bits, ctx->ord);
   N = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

   packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   /* pack exponent vector */
   mpoly_set_monomial(packed_exp, exp, poly->bits, ctx->n, deg, rev);

   /* work out at what index term is */
   exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, maskhi, masklo);

   if (!exists) /* term with that exponent doesn't exist */
      c = 0;
   else  /* term with that monomial exists */
      c = fmpz_mpoly_get_coeff_si(poly, index, ctx);

   TMP_END; 

   return c;
}
