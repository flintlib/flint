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

#define ALLOC_PER_VAR ((FLINT_BITS+4)/3)

char *
_fmpz_mpoly_get_str_pretty(const fmpz * coeffs, const ulong * exps, slong len,
                        const char ** x_in, slong bits, const mpoly_ctx_t mctx)
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

    bound = 1;
    for (i = 0; i < len; i++)
        bound += fmpz_sizeinbase(coeffs + i, 10) + 1;

    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);
    mpoly_degrees_ffmpz((fmpz *) exponents, exps, len, bits, mctx);

    for (i = 0; i < mctx->nvars; i++)
        bound += (fmpz_sizeinbase(exponents + i, 10) + strlen(x[i]) + 3)*len;

    str = flint_malloc(bound);
    off = 0;

    for (i = 0; i < len; i++)
    {
        if (fmpz_sgn(coeffs + i) > 0 && i != 0)
            str[off++] = '+';
        if (coeffs[i] == -WORD(1))
            str[off++] = '-';
        if (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1))
        {
            if (!COEFF_IS_MPZ(coeffs[i]))
                off += flint_sprintf(str + off, "%wd", coeffs[i]);
            else
                off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(coeffs[i]));
        }

        mpoly_get_monomial_ffmpz(exponents, exps + N*i, bits, mctx);

        first = 1;

        for (j = 0; j < mctx->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));
            if (cmp > 0)
            {
                if (!first || (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1)))
                    off += flint_sprintf(str + off, "*");
                off += flint_sprintf(str + off, "%s^", x[j]);

                if (!COEFF_IS_MPZ(exponents[j]))
                    off += flint_sprintf(str + off, "%wd", exponents[j]);
                else
                    off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(exponents[j]));

                first = 0;
            }
            else if (cmp == 0)
            {
                if (!first || (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1)))
                    off += flint_sprintf(str + off, "*");
                off += flint_sprintf(str + off, "%s", x[j]);
                first = 0;
            }
        }

        if (mpoly_monomial_is_zero(exps + i*N, N) && (coeffs[i] == WORD(1) || coeffs[i] == -WORD(1)))
            off += flint_sprintf(str + off, "1");
    }

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;

    return str;
}

char *
fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)
{
   return _fmpz_mpoly_get_str_pretty(poly->coeffs, poly->exps,
                             poly->length, x, poly->bits, ctx->minfo);
}
