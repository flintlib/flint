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

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#include "fmpq_poly.h"
#include "padic_poly.h"

void padic_poly_get_fmpq_poly(fmpq_poly_t rop, 
                              const padic_poly_t op, const padic_ctx_t ctx)
{
    const long len = op->length;

    if (padic_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, len);
    _fmpq_poly_set_length(rop, len);

    if (op->val == 0)
    {
        _fmpz_vec_set(rop->coeffs, op->coeffs, len);
        fmpz_set_ui(rop->den, 1);
    }
    else if (op->val == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(rop->coeffs, op->coeffs, len, ctx->p);
        fmpz_set_ui(rop->den, 1);
    }
    else if (op->val > 1)
    {
        fmpz_t x;

        fmpz_init(x);
        fmpz_pow_ui(x, ctx->p, op->val);

        _fmpz_vec_scalar_mul_fmpz(rop->coeffs, op->coeffs, len, x);
        fmpz_set_ui(rop->den, 1);

        fmpz_clear(x);
    }
    else
    {
        _fmpz_vec_set(rop->coeffs, op->coeffs, len);
        fmpz_pow_ui(rop->den, ctx->p, - ((ulong) op->val));
    }
}

