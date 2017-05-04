/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

int
_fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly, 
                        const ulong * exps, slong len, const char ** x,
                                slong bits, slong n, int deg, int rev, slong N)
{
   slong i, j;
   ulong * degs;
   int r, first;

   TMP_INIT;

   if (len == 0)
   {
        r = fputc('0', file);
        r = (r != EOF) ? 1 : EOF;
        return r;
   }

   TMP_START;

   degs = (ulong *) TMP_ALLOC((n - deg)*sizeof(ulong));
   
   r = 1;
   for (i = len - 1; r > 0 && i >= 0; i--)
   {
      if (fmpz_sgn(poly + i) > 0 && i != len - 1)
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
      {
         r = fmpz_fprint(file, poly + i);
      }

      if (r > 0)
      {
         /* code below expects monomials in opposite order */
         mpoly_get_monomial(degs, exps + i*N, bits, n, deg, 1 - rev);
      }

      first = 1;

      for (j = 0; r > 0 && j < n - deg; j++)
      {
         if (degs[j] > 1)
         {
            if (!first || (poly[i] != WORD(1) && poly[i] != WORD(-1)))
            {
               r = fputc('*', file);
               r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
               r = flint_fprintf(file, "*%s^%wd", x[j], degs[j]);
            first = 0;
         }
         if (degs[j] == 1)
         {
            if (!first || (poly[i] != WORD(1) && poly[i] != WORD(-1)))
            {
               r = fputc('*', file);
               r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
               r = flint_fprintf(file, "*%s", x[j]);
            first = 0;
         }
      }

      if (r > 0 && mpoly_monomial_is_zero(exps + i*N, N) &&
                  (poly[i] == WORD(1) || poly[i] == WORD(-1)))
      {
         r = flint_fprintf(file, "1");
      } 
   }     
   
   TMP_END;

   return r;
}

int
fmpz_mpoly_fprint_pretty(FILE * file, const fmpz_mpoly_t poly,
                                   const char ** x, const fmpz_mpoly_ctx_t ctx)
{
   int deg, rev;

   slong N = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

   degrev_from_ord(deg, rev, ctx->ord);

   return _fmpz_mpoly_fprint_pretty(file, poly->coeffs, poly->exps,
                             poly->length, x, poly->bits, ctx->n, deg, rev, N);

   return 0; 
}
