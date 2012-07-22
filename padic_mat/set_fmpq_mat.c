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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <limits.h>

#include "padic_mat.h"

void padic_mat_set_fmpq_mat(padic_mat_t B, 
                            const fmpq_mat_t A, const padic_ctx_t ctx)
{
    if (!(A->r == padic_mat_nrows(B) && A->c == padic_mat_ncols(B)))
    {
        printf("ERROR (padic_mat_set_fmpq_mat).  Incompatible dimensions.\n");
        abort();
    }

    if (!fmpq_mat_is_empty(A))
    {
        long i, j, m = LONG_MAX, v, w;
        fmpz_t f, g, s, t;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(s);
        fmpz_init(t);

        /* Find min valuation m */
        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                if (!fmpq_is_zero(fmpq_mat_entry(A, i, j)))
                {
                    v = fmpz_remove(t, fmpq_mat_entry_num(A, i, j), ctx->p);
                    w = fmpz_remove(t, fmpq_mat_entry_den(A, i, j), ctx->p);
                    m = FLINT_MIN(m, v - w);
                }

        if (m >= ctx->N)
        {
            padic_mat_zero(B);
        }
        else
        {
            for (i = 0; i < A->r; i++)
                for (j = 0; j < A->c; j++)
                    if (fmpq_is_zero(fmpq_mat_entry(A, i, j)))
                    {
                        fmpz_zero(padic_mat_unit(B, i, j));
                    }
                    else
                    {
                        v = fmpz_remove(s, fmpq_mat_entry_num(A, i, j), ctx->p);
                        w = fmpz_remove(t, fmpq_mat_entry_den(A, i, j), ctx->p);
                        if (v - w >= ctx->N)
                        {
                            fmpz_zero(padic_mat_unit(B, i, j));
                        }
                        else
                        {
                            fmpz_pow_ui(f, ctx->p, (v - w) - m);
                            fmpz_pow_ui(g, ctx->p, ctx->N - (v - w));
                            _padic_inv(t, t, ctx->p, ctx->N - (v - w));
                            fmpz_mul(padic_mat_unit(B, i, j), s, t);
                            fmpz_mod(padic_mat_unit(B, i, j), padic_mat_unit(B, i, j), g);
                            fmpz_mul(padic_mat_unit(B, i, j), padic_mat_unit(B, i, j), f);
                        }
                    }
            padic_mat_val(B) = m;
        }

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(s);
        fmpz_clear(t);
    }
}

