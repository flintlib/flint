/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz_mpoly.h"

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

int
fmpz_mpoly_fprint_pretty(FILE * file, const fmpz_mpoly_t poly,
                                   const char ** x, const fmpz_mpoly_ctx_t ctx)
{
   return _fmpz_mpoly_fprint_pretty(file, poly->coeffs, poly->exps,
                                      poly->length, x, poly->bits, ctx->minfo);
}
