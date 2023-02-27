/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpq_mpoly.h"


int
fmpq_mpoly_fprint_pretty(FILE * file, const fmpq_mpoly_t A,
                               const char ** x_in, const fmpq_mpoly_ctx_t qctx)
{
    int r = 0;
    fmpq_t c;
    slong i, j, N;
    fmpz * exponents;
    const fmpz_mpoly_struct * poly = A->zpoly;
    const mpoly_ctx_struct * mctx = qctx->zctx->minfo;
    char ** x = (char **) x_in;
    TMP_INIT;

    TMP_START;
    fmpq_init(c);

    N = mpoly_words_per_exp(poly->bits, mctx);

    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);

    r = 0;

    if (poly->length == 0)
    {
        r = fputc('0', file);
        goto cleanup;
    }

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(mctx->nvars*sizeof(char *));
        for (i = 0; i < mctx->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(22*sizeof(char));
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }

    for (i = 0; i < poly->length; i++)
    {
        int first = 1;

        fmpq_mul_fmpz(c, A->content, poly->coeffs + i);

        r = flint_fprintf(file, (fmpq_sgn(c) >= 0) ? (i > 0 ? " + " : "")
                                                   : (i > 0 ? " - " : "-") );

        fmpq_abs(c, c);
        if (!fmpq_is_one(c))
        {
            first = 0;
            fmpq_fprint(file, c);
        }

        mpoly_get_monomial_ffmpz(exponents, poly->exps + N*i, poly->bits, mctx);

        for (j = 0; j < mctx->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, UWORD(1));
            if (cmp < 0)
                continue;

            if (!first)
            {
                r = fputc('*', file);
            }

            r = flint_fprintf(file, "%s", x[j]);
            if (cmp > 0)
            {
                r = fputc('^', file);
                r = fmpz_fprint(file, exponents + j);
            }

            first = 0;
        }

        if (first)
        {
            r = flint_fprintf(file, "1");
        }
    }

cleanup:

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    fmpq_clear(c);
    TMP_END;
    return r;
}



