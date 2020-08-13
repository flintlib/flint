/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

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
        flint_printf("ERROR (_padic_mat_fprint).  Mode PADIC_SERIES not implemented yet.\n");
        flint_abort();
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
        flint_printf("ERROR (_padic_mat_fprint).  Unknown print mode.\n");
        flint_abort();
    }

    return 1;
}

