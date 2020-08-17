/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fq_zech_mpoly.h"

int fq_zech_mpoly_fprint_pretty(FILE * file, const fq_zech_mpoly_t A,
                             const char ** x_in, const fq_zech_mpoly_ctx_t ctx)
{
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    int r = 0;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        r = (EOF != fputc('0', file));
        return r;
    }

#define CHECK_r if (r <= 0) goto done;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(ctx->minfo->nvars*sizeof(char *));
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i+1);
        }
    }

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            r = flint_fprintf(file, " + ");
            CHECK_r
        }

        r = flint_fprintf(file, "(");
        CHECK_r
        r = fq_zech_fprint_pretty(file, A->coeffs + i, ctx->fqctx);
        CHECK_r
        r = flint_fprintf(file, ")");
        CHECK_r

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, 1);

            if (cmp > 0)
            {
                r = flint_fprintf(file, "*%s^", x[j]);
                CHECK_r
                r = fmpz_fprint(file, exponents + j);
                CHECK_r
            }
            else if (cmp == 0)
            {
                r = flint_fprintf(file, "*%s", x[j]);
                CHECK_r
            }
        }
    }

done:

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;

    return r;
}
