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

void padic_mat_scalar_div_fmpz(padic_mat_t B, 
                               const padic_mat_t A, const fmpz_t c, 
                               const padic_ctx_t ctx)
{
    if (padic_mat_is_empty(B))
    {
        return;
    }

    if (fmpz_is_zero(c))
    {
        printf("ERROR (padic_mat_scalar_div_fmpz).  c is zero.\n");
        abort();
    }

    if (padic_mat_is_zero(A))
    {
        padic_mat_zero(B);
    }
    else
    {
        fmpz_t d;
        long v;

        fmpz_init(d);
        v = fmpz_remove(d, c, ctx->p);

        if (padic_mat_val(A) - v >= ctx->N)
        {
            padic_mat_zero(B);
        }
        else
        {
            _padic_inv(d, d, ctx->p, ctx->N - padic_mat_val(A) + v);

            fmpz_mat_scalar_mul_fmpz(padic_mat(B), padic_mat(A), d);
            padic_mat_val(B) = padic_mat_val(A) - v;

            _padic_mat_reduce(B, ctx);
        }

        fmpz_clear(d);
    }
}

