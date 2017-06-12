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

void fmpz_mpoly_randtest(fmpz_mpoly_t poly, flint_rand_t state,
    slong length, slong exp_bound, slong coeff_bits, const fmpz_mpoly_ctx_t ctx)
{
   slong i, j, vars;
   fmpz_t c;
   ulong * exp;
   int deg, rev;
   TMP_INIT;

   TMP_START;

   fmpz_init(c);

   degrev_from_ord(deg, rev, ctx->ord);

   vars = ctx->n - deg;

   exp = (ulong *) TMP_ALLOC(vars*sizeof(ulong));

   fmpz_mpoly_zero(poly, ctx);

   for (i = 0; i < length; i++)
   {
      for (j = 0; j < vars; j++)
         exp[j] = n_randint(state, exp_bound);

      fmpz_randtest(c, state, coeff_bits);

      fmpz_mpoly_set_term_fmpz(poly, exp, c, ctx);
   }    

   fmpz_clear(c);
}
