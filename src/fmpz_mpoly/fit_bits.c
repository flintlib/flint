/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_fit_bits(fmpz_mpoly_t A,
                                  flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)
{
   if (A->bits < bits)
   {
      if (A->alloc != 0)
      {
         slong N = mpoly_words_per_exp(bits, ctx->minfo);
         ulong * t = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
         mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, ctx->minfo);
         flint_free(A->exps);
         A->exps = t;
      }

      A->bits = bits;
   }
}

slong fmpz_mpoly_max_bits(const fmpz_mpoly_t A)
{
    return _fmpz_vec_max_bits(A->coeffs, A->length);
}
