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

#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "qadic.h"

/*
    Reduces rop modulo p^e.
 */

void qadic_scalar_mod_ppow(qadic_t rop, const qadic_t op, long e, 
                           const qadic_ctx_t ctx)
{
    e = e - op->val;

    if (e <= 0 || op->length == 0)
    {
        qadic_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, op->length);

        if (e == 1)
        {
            _fmpz_vec_scalar_mod_fmpz(rop->coeffs, 
                                      op->coeffs, op->length, (&ctx->pctx)->p);
        }
        else 
        {
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, e, &ctx->pctx);

            _fmpz_vec_scalar_mod_fmpz(rop->coeffs, 
                                      op->coeffs, op->length, pow);

            if (alloc)
                fmpz_clear(pow);
        }

        _padic_poly_set_length(rop, op->length);
        _padic_poly_normalise(rop);
        rop->val = op->val;
    }
}
