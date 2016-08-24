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

void _exp_get_degrees1(ulong * expvec, ulong v, slong bits,
                                                     slong n, int deg, int rev)
{
   slong i, k = FLINT_BITS/bits;
   ulong mask = bits == FLINT_BITS ? ~UWORD(0) : (UWORD(1) << bits) - UWORD(1);

   v >>= (k - n)*bits;

   if (rev)
   {
      for (i = 0; i < n - deg; i++)
      {
         expvec[i] = (v & mask);
         v >>= bits;
      }
   } else
   {
      for (i = 0; i < n - deg; i++)
      {
         expvec[n - deg - i - 1] = (v & mask);
         v >>= bits;
      }
   }
}

int
_fmpz_mpoly_fprint_pretty1(FILE * file, fmpz * poly, ulong * exps, slong len,
                              char ** x, slong bits, slong n, int deg, int rev)
{
   slong i, j;
   ulong * degs;
   ulong v;
   int r, first;

   if (len == 0)
   {
        r = fputc('0', file);
        r = (r != EOF) ? 1 : EOF;
        return r;
   }

   degs = (ulong *) flint_malloc((n - deg)*sizeof(ulong));
   
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
         v = exps[i];
         _exp_get_degrees1(degs, v, bits, n, deg, rev);
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

      if (r > 0 && v == 0 && (poly[i] == WORD(1) || poly[i] == WORD(-1)))
      {
         r = flint_fprintf(file, "1");
      } 
   }     
   
   flint_free(degs);

   return r;
}

int
fmpz_mpoly_fprint_pretty(FILE * file, fmpz_mpoly_t poly, char ** x, 
                                                          fmpz_mpoly_ctx_t ctx)
{
   int deg, rev;

   if (ctx->N == 1)
   {
      degrev_from_ord(deg, rev, ctx->ord);

      return _fmpz_mpoly_fprint_pretty1(file, poly->coeffs, poly->exps,
                                 poly->length, x, ctx->bits, ctx->n, deg, rev);
   } else
      flint_throw(FLINT_ERROR, "Not implemented yet");

   return 0; 
}
