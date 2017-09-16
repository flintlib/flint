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

void fmpz_mpoly_set_term_si(fmpz_mpoly_t poly,
                        ulong const * exp, slong c, const fmpz_mpoly_ctx_t ctx)
{
   slong i, N, index, exp_bits;
   ulong maskhi, masklo;
   ulong * packed_exp;
   int exists, deg, rev;

   TMP_INIT;

   TMP_START;

   degrev_from_ord(deg, rev, ctx->ord);

   /* compute how many bits are required to represent exp */
   exp_bits = mpoly_exp_bits(exp, ctx->n, deg);
   if (exp_bits > FLINT_BITS)
       flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_set_term_si");

   /* reallocate the number of bits of the exponents of the polynomial */
   exp_bits = mpoly_optimize_bits(exp_bits, ctx->n);
   fmpz_mpoly_fit_bits(poly, exp_bits, ctx);

   masks_from_bits_ord(maskhi, masklo, poly->bits, ctx->ord);
   N = words_per_exp(ctx->n, poly->bits);

   packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   /* pack exponent vector */
   mpoly_set_monomial(packed_exp, exp, poly->bits, ctx->n, deg, rev);
       
   /* work out at what index term should be placed */
   exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, maskhi, masklo);

   if (!exists) /* term with that exponent doesn't exist */
   {
      if (c != 0) /* only set if coeff is nonzero */
      {       
         fmpz_mpoly_fit_length(poly, poly->length + 1, ctx);

         /* shift coeffs and exps by one to make space */
         for (i = poly->length; i >= index + 1; i--)
         {
            fmpz_set(poly->coeffs + i, poly->coeffs + i - 1);
            mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i - 1), N);
         }
      
         /* set exponent */
         mpoly_monomial_set(poly->exps + N*index, packed_exp, N);

         poly->length++; /* safe because length is increasing */

         /* set coeff */
         fmpz_mpoly_set_coeff_si(poly, index, c, ctx);
      }
   } else if (c == 0) /* zero coeff, remove term */
   {
      /* shift coeffs and exps by one to make space */
      for (i = index; i < poly->length - 1; i++)
      {
         fmpz_set(poly->coeffs + i, poly->coeffs + i + 1);
         mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i + 1), N);
      }

      _fmpz_mpoly_set_length(poly, poly->length - 1, ctx);
   } else /* term with that monomial exists, coeff is nonzero */
      fmpz_mpoly_set_coeff_si(poly, index, c, ctx);  

   TMP_END; 
}
