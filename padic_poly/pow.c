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

void _padic_poly_pow(fmpz *rop, long *rval, 
                     const fmpz *op, long val, long len, ulong e,
                     const padic_ctx_t ctx)
{
    fmpz_t pow;
    int alloc;

    *rval = (long) e * val;

    alloc = _padic_ctx_pow_ui(pow, ctx->N - *rval, ctx);

    _fmpz_mod_poly_pow(rop, op, len, e, pow);

    if (alloc)
        fmpz_clear(pow);
}

void padic_poly_pow(padic_poly_t rop, const padic_poly_t op, ulong e, 
                    const padic_ctx_t ctx)
{
    if (e == 0)
    {
        padic_poly_one(rop, ctx);
    }
    else if (op->length == 0 || (long) e * op->val >= ctx->N)
    {
        padic_poly_zero(rop);
    }
    else if (e == 1)
    {
        padic_poly_set(rop, op);
        padic_poly_reduce(rop, ctx);
    }
    else
    {
        const long rlen = (long) e * (op->length - 1) + 1;
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(rlen);
        }
        else
        {
            padic_poly_fit_length(rop, rlen);
            t = rop->coeffs;
        }

        _padic_poly_pow(t, &(rop->val), 
                        op->coeffs, op->val, op->length, e, ctx);

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = rlen;
        }
        _padic_poly_set_length(rop, rlen);
        _padic_poly_normalise(rop);
    }
}

