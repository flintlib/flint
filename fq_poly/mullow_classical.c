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
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_mullow_classical(fq_struct *rop, 
                               const fq_struct *op1, long len1, 
                               const fq_struct *op2, long len2, long n, 
                               const fq_ctx_t ctx)
{
    if ((len1 == 1 && len2 == 1) || n == 1)
    {
        fq_mul(rop, op1, op2, ctx);
    }
    else
    {
        long i;

        _fq_poly_scalar_mul_fq(rop, op1, FLINT_MIN(len1, n), op2, ctx);

        if (n > len1)
            _fq_poly_scalar_mul_fq(rop + len1, 
                                   op2 + 1, n - len1, op1 + len1 - 1, ctx);

        for (i = 0; i < FLINT_MIN(len1, n) - 1; i++)
            _fq_poly_scalar_addmul_fq(rop + i + 1, 
                                      op2 + 1, FLINT_MIN(len2, n - i) - 1, 
                                      op1 + i, ctx);
    }
}

void fq_poly_mullow_classical(fq_poly_t rop, 
    const fq_poly_t op1, const fq_poly_t op2, long n, const fq_ctx_t ctx)
{
    const long len = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0 || n == 0)
    {
        fq_poly_zero(rop);
        return;
    }

    if (n > len)
        n = len;

    if (rop == op1 || rop == op2)
    {
        fq_poly_t t;

        fq_poly_init2(t, n);
        _fq_poly_mullow_classical(t->coeffs, op1->coeffs, op1->length, 
                                             op2->coeffs, op2->length, n, ctx);
        fq_poly_swap(rop, t);
        fq_poly_clear(t);
    }
    else
    {
        fq_poly_fit_length(rop, n);
        _fq_poly_mullow_classical(rop->coeffs, op1->coeffs, op1->length, 
                                               op2->coeffs, op2->length, n, ctx);
    }

    _fq_poly_set_length(rop, n);
    _fq_poly_normalise(rop);
}
