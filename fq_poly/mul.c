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

void _fq_poly_mul(fq_struct *rop, 
                  const fq_struct *op1, long len1, 
                  const fq_struct *op2, long len2, const fq_ctx_t ctx)
{
    const long d = fq_ctx_degree(ctx);

    if (FLINT_MAX(len1, len2) < 6)
    {
        _fq_poly_mul_classical(rop, op1, len1, op2, len2, ctx);
    }
    else if (d < 4)
    {
        _fq_poly_mul_reorder(rop, op1, len1, op2, len2, ctx);
    }
    else
    {
        _fq_poly_mul_KS(rop, op1, len1, op2, len2, ctx);
    }
}

void fq_poly_mul(fq_poly_t rop, 
    const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)
{
    const long len1 = op1->length;
    const long len2 = op2->length;
    const long rlen = op1->length + op2->length - 1;

    if (len1 == 0 || len2 == 0)
    {
        fq_poly_zero(rop);
        return;
    }

    if (rop == op1 || rop == op2)
    {
        fq_poly_t t;

        fq_poly_init2(t, rlen);
        _fq_poly_mul(t->coeffs, op1->coeffs, len1, 
                                op2->coeffs, len2, ctx);
        fq_poly_swap(rop, t);
        fq_poly_clear(t);
    }
    else
    {
        fq_poly_fit_length(rop, rlen);
        _fq_poly_mul(rop->coeffs, op1->coeffs, len1, 
                                  op2->coeffs, len2, ctx);
    }

    _fq_poly_set_length(rop, rlen);
}
