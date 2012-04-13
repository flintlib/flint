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

/*
    Assumes that len > 0.
 */

void _padic_poly_get_fmpq_poly(fmpz *rop, fmpz_t den, 
                               const fmpz *op, long val, long len, 
                               const fmpz_t p)
{
    if (val == 0)
    {
        _fmpz_vec_set(rop, op, len);
        fmpz_one(den);
    }
    else if (val == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(rop, op, len, p);
        fmpz_one(den);
    }
    else if (val > 1)
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_pow_ui(t, p, val);

        _fmpz_vec_scalar_mul_fmpz(rop, op, len, t);
        fmpz_one(den);

        fmpz_clear(t);
    }
    else
    {
        _fmpz_vec_set(rop, op, len);
        fmpz_pow_ui(den, p, -val);
    }
}

void padic_poly_get_fmpq_poly(fmpq_poly_t rop, 
                              const padic_poly_t op, const padic_ctx_t ctx)
{
    const long len = op->length;

    if (len == 0)
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, len);
        _padic_poly_get_fmpq_poly(rop->coeffs, rop->den, 
                                  op->coeffs, op->val, op->length, ctx->p);
        _fmpq_poly_set_length(rop, len);
    }
}

