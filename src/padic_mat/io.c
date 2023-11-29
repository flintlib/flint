/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpq.h"
#include "padic_mat.h"

/* printing *******************************************************************/

int padic_mat_fprint(FILE * file, const padic_mat_t A, const padic_ctx_t ctx)
{
    const slong r = padic_mat(A)->r;
    const slong c = padic_mat(A)->c;

    if (padic_mat_is_empty(A))
    {
        flint_fprintf(file, "%wd %wd\n", r, c);
        return 1;
    }

    if (ctx->mode == PADIC_TERSE)
    {
        slong i, j, v;
        fmpz_t s, t;

        fmpz_init(s);
        fmpz_init(t);

        flint_fprintf(file, "%wd %wd ", r, c);

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
            {
                flint_fprintf(file, " ");

                if (fmpz_is_zero(padic_mat_entry(A, i, j)))
                {
                    flint_fprintf(file, "0");
                }
                else
                {
                    v = A->val
                      + fmpz_remove(t, padic_mat_entry(A, i, j), ctx->p);

                    if (v == 0)
                    {
                        fmpz_fprint(file, t);
                    }
                    else if (v > 0)
                    {
                        fmpz_pow_ui(s, ctx->p, v);
                        fmpz_mul(t, s, t);
                        fmpz_fprint(file, t);
                    }
                    else  /* v < 0 */
                    {
                        fmpz_pow_ui(s, ctx->p, -v);
                        _fmpq_fprint(file, t, s);
                    }
                }
            }

        fmpz_clear(s);
        fmpz_clear(t);
    }
    else if (ctx->mode == PADIC_SERIES)
    {
        flint_throw(FLINT_ERROR, "ERROR (_padic_mat_fprint).  Mode PADIC_SERIES not implemented yet.\n");
    }
    else if (ctx->mode == PADIC_VAL_UNIT)
    {
        slong i, j, v;
        fmpz_t t;

        fmpz_init(t);

        flint_fprintf(file, "%wd %wd ", r, c);

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
            {
                flint_fprintf(file, " ");

                if (fmpz_is_zero(padic_mat_entry(A, i, j)))
                {
                    flint_fprintf(file, "0");
                }
                else
                {
                    v = A->val
                      + fmpz_remove(t, padic_mat_entry(A, i, j), ctx->p);

                    if (v == 0)
                    {
                        fmpz_fprint(file, t);
                    }
                    else if (v == 1)
                    {
                        fmpz_fprint(file, ctx->p);
                        flint_fprintf(file, "*");
                        fmpz_fprint(file, t);
                    }
                    else
                    {
                        fmpz_fprint(file, ctx->p);
                        flint_fprintf(file, "^%wd*", v);
                        fmpz_fprint(file, t);
                    }
                }
            }

        fmpz_clear(t);
    }
    else
    {
        flint_throw(FLINT_ERROR, "ERROR (_padic_mat_fprint).  Unknown print mode.\n");
    }

    return 1;
}

int padic_mat_fprint_pretty(FILE * file, const padic_mat_t A, const padic_ctx_t ctx)
{
    const slong r = padic_mat(A)->r;
    const slong c = padic_mat(A)->c;

    slong i, j, v;
    fmpz_t u;

    fmpz_init(u);

    fputc('[', file);
    for (i = 0; i < r; i++)
    {
        fputc('[', file);
        for (j = 0; j < c; j++)
        {
            v = A->val + fmpz_remove(u, padic_mat_entry(A, i, j), ctx->p);

            _padic_fprint(file, u, v, ctx);

            if (j != c - 1)
                fputc(' ', file);
        }
        fputc(']', file);
        if (i != r - 1)
            fputc('\n', file);
    }
    fputc(']', file);

    fmpz_clear(u);

    return 1;
}

int padic_mat_print(const padic_mat_t A, const padic_ctx_t ctx) { return padic_mat_fprint(stdout, A, ctx); }
int padic_mat_print_pretty(const padic_mat_t A, const padic_ctx_t ctx) { return padic_mat_fprint_pretty(stdout, A, ctx); }
