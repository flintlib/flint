/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "nmod_mpoly.h"

#define ALLOC_PER_VAR ((FLINT_BITS+4)/3)

static char *
_nmod_mpoly_get_str_pretty(const mp_limb_t * coeff, const ulong * exp, slong len,
           const char ** x_in, slong bits, const mpoly_ctx_t mctx, nmod_t fctx)
{
    char * str, ** x = (char **) x_in, *xtmp;
    slong i, j, N, bound, off;
    fmpz * exponents;
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
        xtmp = (char *) TMP_ALLOC(mctx->nvars * ALLOC_PER_VAR * sizeof(char));
        x = (char **) TMP_ALLOC(mctx->nvars*sizeof(char *));
        for (i = 0; i < mctx->nvars; i++)
        {
            x[i] = xtmp + i * ALLOC_PER_VAR;
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }

    bound = 1 + len * ((FLINT_BIT_COUNT(fctx.n) + 3)/3);

    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(ulong));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);
    mpoly_degrees_ffmpz((fmpz *) exponents, exp, len, bits, mctx);

    for (i = 0; i < mctx->nvars; i++)
        bound += (fmpz_sizeinbase(exponents + i, 10) + strlen(x[i]) + 3)*len;

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
            off += flint_sprintf(str + off, "%wu", coeff[i]);
        }

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, mctx);

        for (j = 0; j < mctx->nvars; j++)
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

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;
    return str;
}

char *
nmod_mpoly_get_str_pretty(const nmod_mpoly_t A, const char ** x,
                                                    const nmod_mpoly_ctx_t ctx)
{
   return _nmod_mpoly_get_str_pretty(A->coeffs, A->exps, A->length,
                                             x, A->bits, ctx->minfo, ctx->mod);
}
