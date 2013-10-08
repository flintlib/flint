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

#include "ulong_extras.h"
#include "fq.h"

/*
    Sets (rop, 2d-1) to the image of (op, len) under the Frobenius operator 
    raised to the e-th power, assuming that neither op nor e are zero.
 */

void _fq_frobenius(fmpz *rop, const fmpz *op, long len, long e, 
                   const fmpz *a, const long *j, long lena, 
                   const fmpz_t p)
{
    const long d = j[lena - 1];

    if (len == 1)  /* op is in Fp, not just Fq */
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, (2*d - 1)  - len);
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_pow_ui(t, p, e);
        _fq_pow(rop, op, len, t, a, j, lena, p);
        fmpz_clear(t);
    }
}

void fq_frobenius(fq_t rop, const fq_t op, long e, const fq_ctx_t ctx)
{
    const long d = fq_ctx_degree(ctx);

    e = e % d;
    if (e < 0)
        e += d;

    if (fq_is_zero(op))
    {
        fq_zero(rop);
    }
    else if (e == 0)
    {
        fq_set(rop, op);
    }
    else
    {
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(2 * d - 1);
        }
        else
        {
            fmpz_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        _fq_frobenius(t, op->coeffs, op->length, e, 
                      ctx->a, ctx->j, ctx->len, fq_ctx_prime(ctx));

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            _fmpz_poly_set_length(rop, d);
        }
        _fmpz_poly_normalise(rop);
    }
}

