/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gr_mpoly.h"

static char * _gr_mpoly_default_vars[8] = {
    "x", "y", "z", "s", "t", "u", "v", "w"
};

/* todo: error handling */
int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A,
                             const char ** x_in, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        gr_stream_write(out, "0");
        return GR_SUCCESS;
    }

    N = mpoly_words_per_exp(bits, mctx);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(mctx->nvars*sizeof(char *));

        if (mctx->nvars <= 8)
        {
            for (i = 0; i < mctx->nvars; i++)
                x[i] = _gr_mpoly_default_vars[i];
        }
        else
        {
            slong per_var = 22;
            char * tmp = TMP_ALLOC(mctx->nvars * per_var);

            for (i = 0; i < mctx->nvars; i++)
            {
                x[i] = tmp + per_var * i;
                flint_sprintf(x[i], "x%wd", i+1);
            }
        }
    }

    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(ulong));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            gr_stream_write(out, " + ");
        }

        gr_stream_write(out, "(");
        gr_write(out, GR_ENTRY(A->coeffs, i, cctx->sizeof_elem), cctx);
        gr_stream_write(out, ")");

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, mctx);

        for (j = 0; j < mctx->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, 1);

            if (cmp > 0)
            {
                gr_stream_write(out, "*");
                gr_stream_write(out, x[j]);
                gr_stream_write(out, "^");
                gr_stream_write_fmpz(out, exponents + j);
            }
            else if (cmp == 0)
            {
                gr_stream_write(out, "*");
                gr_stream_write(out, x[j]);
            }
        }
    }

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;

    return GR_SUCCESS;
}

int gr_mpoly_print_pretty(const gr_mpoly_t A,
                             const char ** x_in, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_mpoly_write_pretty(out, A, x_in, mctx, cctx);
}
