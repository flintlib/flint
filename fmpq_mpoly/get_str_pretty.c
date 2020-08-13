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

#define ALLOC_PER_VAR ((FLINT_BITS+4)/3)

char * fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t A,
                                  const char ** x_in, const fmpq_mpoly_ctx_t qctx)
{
    fmpq_t c;
    slong i, j, N, bound, off;
    fmpz * exponents;
    const fmpz_mpoly_struct * poly = A->zpoly;
    const mpoly_ctx_struct * mctx = qctx->zctx->minfo;
    char ** x = (char **) x_in, *xtmp;
    slong len = A->zpoly->length;
    flint_bitcnt_t bits = A->zpoly->bits;
    char * str;
    TMP_INIT;

    if (poly->length == 0)
    {
        str = (char *) flint_malloc(2);
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    N = mpoly_words_per_exp(poly->bits, mctx);
    TMP_START;
    fmpq_init(c);

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
    {
        fmpq_mul_fmpz(c, A->content, poly->coeffs + i);
        bound += fmpz_sizeinbase(fmpq_numref(c), 10);
        bound += fmpz_sizeinbase(fmpq_denref(c), 10);
        bound += 4;

    }
    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(ulong));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);
    mpoly_degrees_ffmpz((fmpz *) exponents, poly->exps, len, bits, mctx);

    for (i = 0; i < mctx->nvars; i++)
        bound += (fmpz_sizeinbase(exponents + i, 10) + strlen(x[i]) + 3)*len;

    str = flint_malloc(bound);
    off = 0;

    for (i = 0; i < len; i++)
    {
        int first = 1;

        fmpq_mul_fmpz(c, A->content, poly->coeffs + i);

        off += flint_sprintf(str + off, "%s", (fmpq_sgn(c) >= 0)
                                                   ? (i > 0 ? " + " : "")
                                                   : (i > 0 ? " - " : "-") );
        fmpq_abs(c, c);
        if (!fmpq_is_one(c))
        {
            first = 0;
            fmpq_get_str(str + off, 10, c);
            while (str[off])
                off++;
        }

        mpoly_get_monomial_ffmpz(exponents, poly->exps + N*i, bits, mctx);

        for (j = 0; j < mctx->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, UWORD(1));
            if (cmp < 0)
                continue;

            if (!first)
                str[off++] = '*';

            off += flint_sprintf(str + off, "%s", x[j]);
            if (cmp > 0)
            {
                str[off++] = '^';
                if (!COEFF_IS_MPZ(exponents[j]))
                    off += flint_sprintf(str + off, "%wd", exponents[j]);
                else
                    off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(exponents[j]));
            }

            first = 0;
        }

        if (first)
            str[off++] = '1';
    }

    fmpq_clear(c);
    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;

    str[off] = '\0';
    return str;
}
