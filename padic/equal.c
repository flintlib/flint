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

#include "padic.h"

int _padic_equal(const padic_t op1, const padic_t op2)
{
    return (padic_val(op1) == padic_val(op2)) && 
           (fmpz_equal(padic_unit(op1), padic_unit(op2)));
}

int padic_equal(const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
{
    /* Exact equality? */
    if (_padic_equal(op1, op2))
    {
<<<<<<< HEAD:padic/equal.c
        return 1;
=======
        printf("Exception (fmpz_poly_mat_transpose). Incompatible dimensions.\n");
        abort();
>>>>>>> eb2b23e... Updates exceptions in fmpz_poly_factor, fmpz_poly_q, fmpz_vec, and nmod_mat.:fmpz_poly_mat/transpose.c
    }

    /* Cases where either op1 or op2 is zero mod p^N */
    if (padic_is_zero(op1, ctx))
        return padic_is_zero(op2, ctx);
    if (padic_is_zero(op2, ctx))
        return 0;

    if (padic_val(op1) == padic_val(op2))
    {
        fmpz_t d, pow;
        int alloc, ans;

        alloc = _padic_ctx_pow_ui(pow, ctx->N - padic_val(op1), ctx);

        fmpz_init(d);
        fmpz_sub(d, padic_unit(op1), padic_unit(op2));
        ans = fmpz_divisible(d, pow);
        fmpz_clear(d);

        if (alloc)
            fmpz_clear(pow);

        return ans;
    }
    else
    {
        return 0;
    }
}
