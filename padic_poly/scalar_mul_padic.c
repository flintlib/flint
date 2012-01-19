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

#include "padic_poly.h"

void _padic_poly_scalar_mul_padic(fmpz *rop, long *rval, 
                                  const fmpz *op, long val, long len, 
                                  const padic_t c, const padic_ctx_t ctx)
{
    if (_padic_is_zero(c) || val + padic_val(c) >= ctx->N)
    {
        _fmpz_vec_zero(rop, len);
        *rval = 0;
    }
    else
    {
        fmpz_t pow;
        int alloc;

        *rval = val + padic_val(c);

        alloc = _padic_ctx_pow_ui(pow, ctx->N - *rval, ctx);

        _fmpz_vec_scalar_mul_fmpz(rop, op, len, padic_unit(c));

        _fmpz_vec_scalar_mod_fmpz(rop, rop, len, pow);

        if (alloc)
            fmpz_clear(pow);    
    }
}

void padic_poly_scalar_mul_padic(padic_poly_t rop, const padic_poly_t op, 
                                 const padic_t c, const padic_ctx_t ctx)
{
    if (padic_poly_is_zero(op) || _padic_is_zero(c) ||
        op->val + padic_val(c) >= ctx->N)
    {
        padic_poly_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, op->length);
        _padic_poly_set_length(rop, op->length);

        _padic_poly_scalar_mul_padic(rop->coeffs, &(rop->val), 
                                     op->coeffs, op->val, op->length, c, ctx);
    }
}

