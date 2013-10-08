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

void _fq_poly_compose_horner(fq_struct *rop, const fq_struct *op1, long len1, 
                                             const fq_struct *op2, long len2, 
                                             const fq_ctx_t ctx)
{
    if (len1 == 1)
    {
        fq_set(rop, op1 + 0);
    }
    else
    {
        const long alloc = (len1 - 1) * (len2 - 1) + 1;

        long i = len1 - 1, lenr;
        fq_struct *t = _fq_poly_init(alloc);
        
        /*
           Perform the first two steps as one, 
             "res = a(m) * poly2 + a(m-1)".
         */
        {
            lenr = len2;
            _fq_poly_scalar_mul_fq(rop, op2, len2, op1 + i, ctx);
            i--;
            fq_add(rop + 0, rop + 0, op1 + i, ctx);
        }
        while (i--)
        {
            _fq_poly_mul(t, rop, lenr, op2, len2, ctx);
            lenr += len2 - 1;
            _fq_poly_add(rop, t, lenr, op1 + i, 1, ctx);
        }
        
        _fmpz_vec_clear(t, alloc);
    }
}

void fq_poly_compose_horner(fq_poly_t rop, 
                            const fq_poly_t op1, const fq_poly_t op2, 
                            const fq_ctx_t ctx)
{
    const long len1 = op1->length;
    const long len2 = op2->length;
    const long lenr = (len1 - 1) * (len2 - 1) + 1;
    
    if (len1 == 0)
    {
        fq_poly_zero(rop);
    }
    else if (len1 == 1 || len2 == 0)
    {
        fq_poly_set_fq(rop, op1->coeffs + 0);
    }
    else if (rop != op1 && rop != op2)
    {
        fq_poly_fit_length(rop, lenr);
        _fq_poly_compose_horner(rop->coeffs, op1->coeffs, len1, 
                                             op2->coeffs, len2, ctx);
        _fq_poly_set_length(rop, lenr);
        _fq_poly_normalise(rop);
    }
    else
    {
        fq_poly_t t;

        fq_poly_init2(t, lenr);
        _fq_poly_compose_horner(t->coeffs, op1->coeffs, len1,
                                           op2->coeffs, len2, ctx);
        _fq_poly_set_length(t, lenr);
        _fq_poly_normalise(t);
        fq_poly_swap(rop, t);
        fq_poly_clear(t);
    }
}

