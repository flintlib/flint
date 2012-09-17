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

void _fq_poly_sqr_classical(fq_struct *rop, 
                            const fq_struct *op, long len, const fq_ctx_t ctx)
{
    if (len == 1)
    {
        fq_mul(rop, op, op, ctx);
    }
    else
    {
        long i;
        fq_t t;

        fq_init(t);

        _fq_poly_scalar_mul_fq(rop, op, len, op, ctx);

        _fq_poly_scalar_mul_fq(rop + len, op + 1, len - 1, op + len - 1, ctx);

        for (i = 1; i < len - 1; i++)
            _fq_poly_scalar_addmul_fq(rop + i + 1, op + 1, i - 1, op + i, ctx);

        for (i = 1; i < 2 * len - 2; i++)
            fq_add(rop + i, rop + i, rop + i, ctx);

        for (i = 1; i < len - 1; i++)
        {
            fq_sqr(t, op + i, ctx);
            fq_add(rop + 2 * i, rop + 2 * i, t, ctx);
        }
        fq_clear(t);
    }
}

void fq_poly_sqr_classical(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
{
    const long len = 2 * op->length - 1;

    if (op->length == 0)
    {
        fq_poly_zero(rop);
        return;
    }

    if (rop == op)
    {
        fq_poly_t t;

        fq_poly_init2(t, len);
        _fq_poly_sqr_classical(t->coeffs, op->coeffs, op->length, ctx);
        fq_poly_swap(rop, t);
        fq_poly_clear(t);
    }
    else
    {
        fq_poly_fit_length(rop, len);
        _fq_poly_sqr_classical(rop->coeffs, op->coeffs, op->length, ctx);
    }

    _fq_poly_set_length(rop, len);
}
