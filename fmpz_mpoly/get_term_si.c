/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2017 Daniel Schultz

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
   slong c;
   slong N, index, exp_bits;
   ulong * cmpmask;
   ulong * packed_exp;
   int exists;
   TMP_INIT;

   exp_bits = mpoly_exp_bits_required_ui(exp, ctx->minfo);
   if (exp_bits > FLINT_BITS)
       flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_get_term_si");

   if (exp_bits > poly->bits) /* exponent too large to be poly exponent */
       return 0;

   TMP_START;
   
    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);

   packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   /* pack exponent vector */
   mpoly_set_monomial_ui(packed_exp, exp, poly->bits, ctx->minfo);

   /* work out at what index term is */
   exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, cmpmask);

   if (!exists) /* term with that exponent doesn't exist */
      c = 0;
   else  /* term with that monomial exists */
      c = fmpz_mpoly_get_coeff_si(poly, index, ctx);

   TMP_END; 

   return c;
}
