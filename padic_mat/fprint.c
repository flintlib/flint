/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "fmpz_mat.h"
#include "padic_mat.h"

int padic_mat_fprint(FILE * file, const padic_mat_t A, const padic_ctx_t ctx)
{
    const long r = padic_mat(A)->r;
    const long c = padic_mat(A)->c;

    if (padic_mat_is_empty(A))
    {
        fprintf(file, "%ld %ld\n", r, c);
        return 1;
    }

    if (ctx->mode == PADIC_TERSE)
    {
        long i, j, v;
        fmpz_t s, t;

        fmpz_init(s);
        fmpz_init(t);

        fprintf(file, "%ld %ld ", r, c);

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
            {
                fprintf(file, " ");

                if (fmpz_is_zero(padic_mat_unit(A, i, j)))
                {
                    fprintf(file, "0");
                }
                else
                {
                    v = A->val 
                      + fmpz_remove(t, padic_mat_unit(A, i, j), ctx->p);

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
        printf("ERROR (_padic_mat_fprint).  Not implemented yet.\n");
        abort();
    }
    else if (ctx->mode == PADIC_VAL_UNIT)
    {
        long i, j, v;
        fmpz_t t;

        fmpz_init(t);

        fprintf(file, "%ld %ld ", r, c);

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
            {
                fprintf(file, " ");

                if (fmpz_is_zero(padic_mat_unit(A, i, j)))
                {
                    fprintf(file, "0");
                }
                else
                {
                    v = A->val 
                      + fmpz_remove(t, padic_mat_unit(A, i, j), ctx->p);

                    if (v == 0)
                    {
                        fmpz_fprint(file, t);
                    }
                    else if (v == 1)
                    {
                        fmpz_fprint(file, ctx->p);
                        fprintf(file, "*");
                        fmpz_fprint(file, t);
                    }
                    else
                    {
                        fmpz_fprint(file, ctx->p);
                        fprintf(file, "^%ld*", v);
                        fmpz_fprint(file, t);
                    }
                }
            }

        fmpz_clear(t);
    }
    else
    {
        printf("ERROR (_padic_mat_fprint).  Unknown print mode.\n");
        abort();
    }

    return 1;
}

