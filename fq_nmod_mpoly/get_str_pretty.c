/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fq_nmod_mpoly.h"

char *
_fq_nmod_mpoly_get_str_pretty(const fq_nmod_struct * coeff, const ulong * exp,
                                  slong len, const char ** x_in, slong bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    char * str, ** x = (char **) x_in;
    slong i, j, N, bound, off;
    fmpz * exponents;
    int first;
    char ** coeffstrs;
    TMP_INIT;

    if (len == 0)
    {
        str = flint_malloc(2);
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(ctx->minfo->nvars*sizeof(char *));
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }

    coeffstrs = (char **)flint_malloc(len * sizeof(char *));
    
    bound = 1;
    for (i = 0; i < len; i++)
    {
        coeffstrs[i] = fq_nmod_get_str_pretty(coeff + i, ctx->fqctx);
        bound += 5 + strlen(coeffstrs[i]);
    }

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
    mpoly_degrees_ffmpz(exponents, exp, len, bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nvars; i++)
        bound += (2 + fmpz_sizeinbase(exponents + i, 10) + strlen(x[i]) + 3)*len;

    str = flint_malloc(bound);
    off = 0;
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            str[off++] = ' ';
            str[off++] = '+';
            str[off++] = ' ';
        }

        first = fq_nmod_is_one(coeff + i, ctx->fqctx);
        if (!first)
        {
            off += flint_sprintf(str + off, "(%s)", coeffstrs[i]);
        }

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            if (fmpz_is_zero(exponents + j))
                continue;

            if (!first)
            {
                str[off++] = '*';
            }

            if (fmpz_cmp_ui(exponents + j, UWORD(1)) > 0)
            {
                off += flint_sprintf(str + off, "%s^", x[j]);
                if (!COEFF_IS_MPZ(exponents[j]))
                    off += flint_sprintf(str + off, "%wd", exponents[j]);
                else
                    off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(exponents[j]));

            }
            else
            {
                off += flint_sprintf(str + off, "%s", x[j]);
            }
            
            first = 0;
        }

        if (first)
        {
            off += flint_sprintf(str + off, "1");
        }

        FLINT_ASSERT(off < bound);
    }

    for (i = 0; i < len; i++)
    {
        flint_free(coeffstrs[i]);
    }
    flint_free(coeffstrs);

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_clear(exponents + i);
    }

    TMP_END;
    return str;
}

char *
fq_nmod_mpoly_get_str_pretty(const fq_nmod_mpoly_t A, const char ** x,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
   return _fq_nmod_mpoly_get_str_pretty(A->coeffs, A->exps, A->length,
                                                              x, A->bits, ctx);
}
