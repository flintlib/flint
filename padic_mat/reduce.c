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

void _padic_mat_reduce(padic_mat_t mat, const padic_ctx_t ctx)
{
    if (!padic_mat_is_empty(mat) && !padic_mat_is_zero(mat))
    {
        if (mat->val >= ctx->N)
        {
            padic_mat_zero(mat);
        }
        else
        {
            long i;
            fmpz_t pow;

            fmpz_init(pow);
            fmpz_pow_ui(pow, ctx->p, ctx->N - mat->val);
            for (i = 0; i < padic_mat(mat)->r * padic_mat(mat)->c; i++)
            {
                fmpz_mod(padic_mat(mat)->entries + i, 
                         padic_mat(mat)->entries + i, pow);
            }
            fmpz_clear(pow);

            if (padic_mat_is_zero(mat))
            {
                mat->val = 0;
            }
        }
    }
}

void padic_mat_reduce(padic_mat_t mat, const padic_ctx_t ctx)
{
    _padic_mat_canonicalise(mat, ctx);
    _padic_mat_reduce(mat, ctx);
}

