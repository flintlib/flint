/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz_mod_mpoly.h"
#include "fmpz_mpoly.h"

/* printing *******************************************************************/

int fmpz_mod_mpoly_fprint_pretty(FILE * file, const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx) { return _fmpz_mpoly_fprint_pretty(file, A->coeffs, A->exps, A->length, x, A->bits, ctx->minfo); }
int fmpz_mod_mpoly_print_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx) { return fmpz_mod_mpoly_fprint_pretty(stdout, A, x, ctx); }

/* debugging ******************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that no monomial of r is divisible by lm(g)
*/
void fmpz_mod_mpoly_remainder_strongtest(const fmpz_mod_mpoly_t r,
                      const fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_ctx_t ctx)
{
   slong i, N, bits;
   ulong mask = 0;
   ulong * rexp, * gexp;

   bits = FLINT_MAX(r->bits, g->bits);
   N = mpoly_words_per_exp(bits, ctx->minfo);

    if (g->length == 0)
        flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

    if (r->length == 0)
        return;

    rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
    gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
    mpoly_repack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->minfo);
    mpoly_repack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->minfo);

    if (bits <= FLINT_BITS)
        mask = mpoly_overflow_mask_sp(bits);
    else
        mask = 0;

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides)
        {
            flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_remainder_strongtest FAILED i = %wd\n"
                                     "rem %s\n\n"
                                     "den %s\n\n",
                                     i,
                                     fmpz_mod_mpoly_get_str_pretty(r, NULL, ctx),
                                     fmpz_mod_mpoly_get_str_pretty(g, NULL, ctx));
        }
    }

   flint_free(rexp);
   flint_free(gexp);
}
