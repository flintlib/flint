/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz_mpoly.h"

/* printing *******************************************************************/

int
_fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly,
                        const ulong * exps, slong len, const char ** x_in,
                                      flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
   slong i, j, N;
   fmpz * exponents;
   int r, first;
   char ** x = (char **) x_in;
   TMP_INIT;

   if (len == 0)
   {
        r = fputc('0', file);
        r = (r != EOF) ? 1 : EOF;
        return r;
   }

    N = mpoly_words_per_exp(bits, mctx);

   TMP_START;

   if (x == NULL)
   {
      x = (char **) TMP_ALLOC(mctx->nvars*sizeof(char *));

      for (i = 0; i < mctx->nvars; i++)
      {
         x[i] = (char *) TMP_ALLOC(22*sizeof(char));
         flint_sprintf(x[i], "x%wd", i + 1);
      }
   }

    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);

   r = 1;
   for (i = 0; r > 0 && i < len; i++)
   {
      if (fmpz_sgn(poly + i) >= 0 && i != 0)
      {
         r = fputc('+', file);
         r = (r != EOF) ? 1 : EOF;
      }
      if (poly[i] == WORD(-1))
      {
         r = fputc('-', file);
         r = (r != EOF) ? 1 : EOF;
      }
      if (r > 0 && poly[i] != WORD(1) && poly[i] != WORD(-1))
         r = fmpz_fprint(file, poly + i);

      if (r > 0)
         mpoly_get_monomial_ffmpz(exponents, exps + N*i, bits, mctx);

      first = 1;

      for (j = 0; r > 0 && j < mctx->nvars; j++)
      {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));
         if (cmp > 0)
         {
            if (!first || (poly[i] != WORD(1) && poly[i] != WORD(-1)))
            {
               r = fputc('*', file);
               r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
               r = flint_fprintf(file, "%s^", x[j]);
            if (r > 0)
                r = fmpz_fprint(file, exponents + j);

            first = 0;
         }
         else if (cmp == 0)
         {
            if (!first || (poly[i] != WORD(1) && poly[i] != WORD(-1)))
            {
               r = fputc('*', file);
               r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
               r = flint_fprintf(file, "%s", x[j]);
            first = 0;
         }
      }

      if (r > 0 && mpoly_monomial_is_zero(exps + N*i, N) &&
                  (poly[i] == WORD(1) || poly[i] == WORD(-1)))
      {
         r = flint_fprintf(file, "1");
      }
   }

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

   TMP_END;

   return r;
}

int fmpz_mpoly_fprint_pretty(FILE * file, const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx) { return _fmpz_mpoly_fprint_pretty(file, poly->coeffs, poly->exps, poly->length, x, poly->bits, ctx->minfo); }
int _fmpz_mpoly_print_pretty(const fmpz * poly, const ulong * exps, slong len, const char ** x, slong bits, const mpoly_ctx_t mctx) { return _fmpz_mpoly_fprint_pretty(stdout, poly, exps, len, x, bits, mctx); }
int fmpz_mpoly_print_pretty(const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx) { return fmpz_mpoly_fprint_pretty(stdout, A, x, ctx); }

/* debugging ******************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that if c*x^a is a term of r and x^a is divisible by the leading
   monomial of g, then |c| < |leading coefficient of g|
*/
void fmpz_mpoly_remainder_test(const fmpz_mpoly_t r, const fmpz_mpoly_t g,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, N, bits;
   ulong mask = 0;
   ulong * rexp, * gexp;

   bits = FLINT_MAX(r->bits, g->bits);
   N = mpoly_words_per_exp(bits, ctx->minfo);

   if (g->length == 0 )
      flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

   if (r->length == 0 )
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

        if (divides && fmpz_cmpabs(g->coeffs + 0, r->coeffs + i) <= 0)
        {
            flint_throw(FLINT_ERROR, "fmpz_mpoly_remainder_test FAILED i = %wd\n"
                                     "rem %s\n\n"
                                     "den %s\n\n",
                                     i,
                                     fmpz_mpoly_get_str_pretty(r, NULL, ctx),
                                     fmpz_mpoly_get_str_pretty(g, NULL, ctx));
        }
    }

   flint_free(rexp);
   flint_free(gexp);
}


/*
   test that r is a valid remainder upon division by g over Q
   this means that no term of r is divisible by lt(g)
*/
void fmpz_mpoly_remainder_strongtest(const fmpz_mpoly_t r, const fmpz_mpoly_t g,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, bits;
    ulong mask = 0;
    ulong * rexp, * gexp;

    bits = FLINT_MAX(r->bits, g->bits);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    if (g->length == 0 )
        flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

    if (r->length == 0 )
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
            flint_throw(FLINT_ERROR, "fmpz_mpoly_remainder_strongtest FAILED i = %wd\n"
                                     "rem %s\n\n"
                                     "den %s\n\n",
                                     i,
                                     fmpz_mpoly_get_str_pretty(r, NULL, ctx),
                                     fmpz_mpoly_get_str_pretty(g, NULL, ctx));
        }
    }

    flint_free(rexp);
    flint_free(gexp);
}
