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

void _fq_poly_sqr(fq_struct *rop, 
                  const fq_struct *op, long len, const fq_ctx_t ctx)
{
    const long d = fq_ctx_degree(ctx);

    if (len < 6)
    {
        _fq_poly_sqr_classical(rop, op, len, ctx);
    }
    else if (d < 4)
    {
        _fq_poly_sqr_reorder(rop, op, len, ctx);
    }
    else
    {
        _fq_poly_sqr_KS(rop, op, len, ctx);
    }
}

void fq_poly_sqr(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
{
    const long rlen = 2 * op->length - 1;

    if (op->length == 0)
    {
        fq_poly_zero(rop);
        return;
    }

    if (rop == op)
    {
        fq_poly_t t;

        fq_poly_init2(t, rlen);
        _fq_poly_sqr(t->coeffs, op->coeffs, op->length, ctx);
        fq_poly_swap(rop, t);
        fq_poly_clear(t);
    }
    else
    {
        fq_poly_fit_length(rop, rlen);
        _fq_poly_sqr(rop->coeffs, op->coeffs, op->length, ctx);
    }

    _fq_poly_set_length(rop, rlen);
}
