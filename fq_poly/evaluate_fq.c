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

void _fq_poly_evaluate_fq(fq_t rop, const fq_struct *op, long len, 
                                    const fq_t a, const fq_ctx_t ctx)
{
    if (len == 0)
    {
        fq_zero(rop);
    }
    else if (len == 1 || fq_is_zero(a))
    {
        fq_set(rop, op + 0);
    }
    else
    {
        long i = len - 1;
        fq_t t;

        fq_init(t);
        fq_set(rop, op + i);
        for (i = len - 2; i >= 0; i--)
        {
            fq_mul(t, rop, a, ctx);
            fq_add(rop, op + i, t, ctx);
        }
        fq_clear(t);
    }
}

void fq_poly_evaluate_fq(fq_t rop, const fq_poly_t f, const fq_t a, 
                                   const fq_ctx_t ctx)
{
    if (rop == a)
    {
        fq_t t;
        fq_init(t);
        _fq_poly_evaluate_fq(t, f->coeffs, f->length, a, ctx);
        fq_swap(rop, t);
        fq_clear(t);
    }
    else
    {
        _fq_poly_evaluate_fq(rop, f->coeffs, f->length, a, ctx);
    }
}

