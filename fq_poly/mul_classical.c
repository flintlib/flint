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

#include "fq_poly.h"

void _fq_poly_mul_classical(fq_struct *rop, 
                            const fq_struct *op1, long len1, 
                            const fq_struct *op2, long len2, 
                            const fq_ctx_t ctx)
{
    if (len1 == 1 && len2 == 1)
    {
        fq_mul(rop, op1, op2, ctx);
    }
    else
    {
        long i;

        /* Set res[i] = poly1[i]*poly2[0] */
        _fq_poly_scalar_mul_fq(rop, op1, len1, op2, ctx);

        /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
        _fq_poly_scalar_mul_fq(rop + len1, op2 + 1, len2 - 1, op1 + len1 - 1, ctx);

        /* out[i+j] += in1[i]*in2[j] */
        for (i = 0; i < len1 - 1; i++)
            _fq_poly_scalar_addmul_fq(rop + i + 1, op2 + 1, len2 - 1, op1 + i, ctx);
    }
}

void fq_poly_mul_classical(fq_poly_t rop, 
    const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)
{
    const long len = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0)
    {
        fq_poly_zero(rop);
        return;
    }

    if (rop == op1 || rop == op2)
    {
        fq_poly_t t;

        fq_poly_init2(t, len);
        _fq_poly_mul_classical(t->coeffs, op1->coeffs, op1->length, 
                                          op2->coeffs, op2->length, ctx);
        fq_poly_swap(rop, t);
        fq_poly_clear(t);
    }
    else
    {
        fq_poly_fit_length(rop, len);
        _fq_poly_mul_classical(rop->coeffs, op1->coeffs, op1->length, 
                                            op2->coeffs, op2->length, ctx);
    }

    _fq_poly_set_length(rop, len);
}
