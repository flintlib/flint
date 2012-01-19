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

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

void _padic_poly_derivative(fmpz *rop, long *rval, 
                            const fmpz *op, long val, long len, 
                            const padic_ctx_t ctx)
{
    fmpz_t pow;
    int alloc;

    _fmpz_poly_derivative(rop, op, len);
    *rval = val;

    _padic_poly_canonicalise(rop, rval, len - 1, ctx->p);

    alloc = _padic_ctx_pow_ui(pow, ctx->N - *rval, ctx);

    _fmpz_vec_scalar_mod_fmpz(rop, rop, len - 1, pow);

    if (alloc)
        fmpz_clear(pow);
}

void padic_poly_derivative(padic_poly_t rop, 
                           const padic_poly_t op, const padic_ctx_t ctx)
{
    const long len = op->length;

    if (len < 2 || op->val >= ctx->N)
    {
        padic_poly_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, len - 1);
        _padic_poly_derivative(rop->coeffs, &(rop->val), 
                               op->coeffs, op->val, len, ctx);
        _padic_poly_set_length(rop, len - 1);
        _padic_poly_normalise(rop);
    }
}

