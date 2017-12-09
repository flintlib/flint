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
#include "nmod_mpoly.h"

char *
_nmod_mpoly_get_str_pretty(const mp_limb_t * coeff, const ulong * exp, slong len,
                             const char ** x_in, slong bits,
                                const mpoly_ctx_t mctx, const nmodf_ctx_t fctx)
{
    char * str, ** x = (char **) x_in;
    slong i, j, N, bound, off;
    ulong * degs;
    int first;

    TMP_INIT;

    if (len == 0)
    {
        str = flint_malloc(2);
        str[0] = '0';
        str[1] = '\0';
        return str;
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

    bound = 1 + len * ((FLINT_BIT_COUNT(fctx->mod.n) + 3)/3);

    degs = (ulong *) TMP_ALLOC(mctx->nvars*sizeof(ulong));
    mpoly_degrees((slong *) degs, exp, len, bits, mctx);

    for (i = 0; i < mctx->nvars; i++)
    {
        bound += ((FLINT_BIT_COUNT(degs[i]) + 3)/3 + strlen(x[i]) + 3)*len;
    }

    str = flint_malloc(bound);
    off = 0;
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            str[off++] = '+';
        }

        first = (coeff[i] == 1);
        if (!first)
        {
            off += flint_sprintf(str + off, "%wd", coeff[i]);
        }

        mpoly_get_monomial_ui(degs, exp + N*i, bits, mctx);

        for (j = 0; j < mctx->nvars; j++)
        {
            if (degs[j] == 0)
                continue;

            if (!first)
            {
                str[off++] = '*';
            }

            if (degs[j] > 1)
            {
                off += flint_sprintf(str + off, "%s^%wd", x[j], degs[j]);
            } else
            {
                off += flint_sprintf(str + off, "%s", x[j]);
            }
            
            first = 0;
        }

        if (first)
        {
            off += flint_sprintf(str + off, "1");
        }
   }

   TMP_END;
   return str;
}

char *
nmod_mpoly_get_str_pretty(const nmod_mpoly_t poly, const char ** x,
                                                    const nmod_mpoly_ctx_t ctx)
{
   return _nmod_mpoly_get_str_pretty(poly->coeffs, poly->exps, poly->length,
                                       x, poly->bits, ctx->minfo, ctx->ffinfo);
}
