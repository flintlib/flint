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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_sqr_KS(fq_struct *rop, const fq_struct *op, long len, 
                                     const fq_ctx_t ctx)
{
    const long in_len = len;
    const long d = fq_ctx_degree(ctx);
    long bits, i;
    fmpz *f, *g;

    FQ_VEC_NORM(op, len);

    if (!len)
    {
        if (2 * in_len - 1 > 0)
            _fq_poly_zero(rop, 2 * in_len - 1);
        return;
    }

    bits = 2 * fmpz_bits(fq_ctx_prime(ctx)) 
           + FLINT_BIT_COUNT(d) + FLINT_BIT_COUNT(len);

    f = _fmpz_vec_init((2 * len - 1) + len);
    g = f + (2 * len - 1);

    for (i = 0; i < len; i++)
    {
        fmpz_poly_bit_pack(g + i, op + i, bits);
    }

    _fmpz_poly_sqr(f, g, len);

    for (i = 0; i < 2 * len - 1; i++)
    {
        fmpz_poly_bit_unpack_unsigned(rop + i, f + i, bits);
        fq_reduce(rop + i, ctx);
    }

    _fq_poly_zero(rop + (2 * len - 1), 2 * (in_len - len));

    _fmpz_vec_clear(f, (2 * len - 1) + len);
}

void fq_poly_sqr_KS(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
{
    const long len = 2 * op->length - 1;

    if (op->length == 0)
    {
        fq_poly_zero(rop);
    }
    else
    {
        fq_poly_fit_length(rop, len);
        _fq_poly_sqr_KS(rop->coeffs, op->coeffs, op->length, ctx);
        _fq_poly_set_length(rop, len);
    }
}
