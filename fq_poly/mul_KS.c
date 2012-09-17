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

void _fq_poly_mul_KS(fq_struct *rop, const fq_struct *op1, long len1, 
                                     const fq_struct *op2, long len2, 
                                     const fq_ctx_t ctx)
{
    const long in1_len = len1, in2_len = len2;
    const long d = fq_ctx_degree(ctx);
    long bits, i;
    fmpz *f, *g, *h;

    FQ_VEC_NORM(op1, len1);
    FQ_VEC_NORM(op2, len2);

    if (!len1 | !len2)
    {
        if (in1_len + in2_len - 1 > 0)
            _fq_poly_zero(rop, in1_len + in2_len - 1);
        return;
    }

    bits = 2 * fmpz_bits(fq_ctx_prime(ctx)) 
           + FLINT_BIT_COUNT(d) + FLINT_BIT_COUNT(FLINT_MIN(len1, len2));

    f = _fmpz_vec_init((len1 + len2 - 1) + (len1) + (len2));
    g = f + (len1 + len2 - 1);
    h = g + len1;

    for (i = 0; i < len1; i++)
    {
        fmpz_poly_bit_pack(g + i, op1 + i, bits);
    }
    for (i = 0; i < len2; i++)
    {
        fmpz_poly_bit_pack(h + i, op2 + i, bits);
    }

    if (len1 >= len2)
        _fmpz_poly_mul(f, g, len1, h, len2);
    else
        _fmpz_poly_mul(f, h, len2, g, len1);

    for (i = 0; i < len1 + len2 - 1; i++)
    {
        fmpz_poly_bit_unpack_unsigned(rop + i, f + i, bits);
        fq_reduce(rop + i, ctx);
    }

    _fq_poly_zero(rop + (len1 + len2 - 1), (in1_len - len1) + (in2_len - len2));

    _fmpz_vec_clear(f, (len1 + len2 - 1) + (len1) + (len2));
}

void fq_poly_mul_KS(fq_poly_t rop, 
                    const fq_poly_t op1, const fq_poly_t op2, 
                    const fq_ctx_t ctx)
{
    const long len1 = op1->length;
    const long len2 = op2->length;
    const long rlen = len1 + len2 - 1;

    if (len1 == 0 || len2 == 0)
    {
        fq_poly_zero(rop);
    }
    else
    {
        fq_poly_fit_length(rop, rlen);
        _fq_poly_mul_KS(rop->coeffs, op1->coeffs, len1, 
                                     op2->coeffs, len2, ctx);
        _fq_poly_set_length(rop, rlen);
    }
}
