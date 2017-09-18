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

void fmpz_mpoly_set_monomial(fmpz_mpoly_t poly, 
                       slong n, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
   slong exp_bits, N;
   int deg, rev;

   degrev_from_ord(deg, rev, ctx->ord);

   if (n > poly->length)
      flint_throw(FLINT_ERROR, "Invalid index in fmpz_mpoly_set_monomial");

   /* compute how many bits are required to represent exp */
   exp_bits = mpoly_exp_bits(exp, ctx->n, deg);
   if (exp_bits > FLINT_BITS)
       flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_set_monomial");

   /* reallocate the number of bits of the exponents of the polynomial */
   exp_bits = mpoly_optimize_bits(exp_bits, ctx->n);
   fmpz_mpoly_fit_bits(poly, exp_bits, ctx);

   N = words_per_exp(ctx->n, poly->bits);

   fmpz_mpoly_fit_length(poly, n + 1, ctx);

   mpoly_set_monomial(poly->exps + n*N, exp, poly->bits, ctx->n, deg, rev);

   if (n == poly->length)
      _fmpz_mpoly_set_length(poly, n + 1, ctx);
}
