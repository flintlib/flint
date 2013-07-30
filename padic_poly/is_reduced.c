/*============================================================================

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

#include "padic_poly.h"

int _padic_poly_is_reduced(const fmpz *op, slong val, slong len, slong N, 
                           const padic_ctx_t ctx)
{
    slong w;

    if (len == 0)
    {
        return (val == 0);
    }
    
    w = _fmpz_vec_ord_p(op, len, ctx->p);

    if (w != 0 || val >= N)
    {
        return 0;
    }

    {
        fmpz_t pow;
        int r, alloc;
        slong i;

        alloc = _padic_ctx_pow_ui(pow, N - val, ctx);

        r = 1;
        for (i = 0; (i < len) && (r); i++)
            if (fmpz_sgn(op + i) < 0 || fmpz_cmp(op + i, pow) >= 0)
                r = 0;

        if (alloc)
            fmpz_clear(pow);

        return r;
    }
}

int padic_poly_is_reduced(const padic_poly_t op, const padic_ctx_t ctx)
{
    return _padic_poly_is_reduced(op->coeffs, op->val, op->length, op->N, ctx);
}

