/*
    Copyright (C) 2017 Daniel Schultz

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
#include "nmod_mpoly.h"

int
_nmod_mpoly_fprint_pretty(FILE * file, const ulong * poly, const ulong * exps,
        slong len, const char ** x_in,  slong bits,
            slong n, int deg, int rev, slong N, const nmodf_ctx_t fctx)
{
    slong i, j, nvars;
    ulong * degs;
    int r = 0, first;
    char ** x = (char **) x_in;

    TMP_INIT;

    if (len == 0)
    {
        r = fputc('0', file);
        r = (r != EOF) ? 1 : EOF;
        return r;
    }

    TMP_START;

    nvars = n - deg;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(nvars*sizeof(char *));
        for (i = 0; i < nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }

    degs = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            r = fputc('+', file);
            r = (r != EOF) ? 1 : EOF;
            if (r <= 0) goto done;
        }

        first = nmodf_is_one(poly + i*fctx->deg, fctx);
        if (!first)
        {
            r = nmodf_fprintf(file, poly + i*fctx->deg, fctx);
            if (r <= 0) goto done;
        }

        mpoly_get_monomial(degs, exps + i*N, bits, n, deg, rev);

        for (j = 0; j < nvars; j++)
        {
            if (degs[j] == 0)
                continue;

            if (!first)
            {
                r = fputc('*', file);
                r = (r != EOF) ? 1 : EOF;
                if (r <= 0) goto done;
            }
            if (degs[j] > 1)
                r = flint_fprintf(file, "%s^%wd", x[j], degs[j]);
            else
                r = flint_fprintf(file, "%s", x[j]);
            if (r <= 0) goto done;
            
            first = 0;
        }

        if (first)
        {
            r = flint_fprintf(file, "1");
            if (r <= 0) goto done;
        }
   }
   
done:
   TMP_END;
   return r;
}

int
nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
   int deg, rev;
   slong N = words_per_exp(ctx->n, poly->bits);
   degrev_from_ord(deg, rev, ctx->ord);

   return _nmod_mpoly_fprint_pretty(file, poly->coeffs, poly->exps,
                poly->length, x, poly->bits, ctx->n, deg, rev, N, ctx->ffinfo);

}
