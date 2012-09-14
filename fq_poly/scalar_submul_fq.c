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

void _fq_poly_scalar_submul_fq(fq_struct *rop, 
    const fq_struct *op, long len, const fq_t x, const fq_ctx_t ctx)
{
    if (fq_is_zero(x))
        return;

    if (fq_is_one(x))
    {
        _fq_poly_sub(rop, rop, len, op, len, ctx);
    }
    else
    {
        long i;
        fq_t t;

        fq_init(t);

        for (i = 0; i < len; i++)
        {
            fq_mul(t, op + i, x, ctx);
            fq_sub(rop + i, rop + i, t, ctx);
        }

        fq_clear(t);
    }
}

void fq_poly_scalar_submul_fq(fq_poly_t rop, 
    const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)
{
    if (!(fq_is_zero(x) || fq_poly_is_zero(op)))
    {
        fq_poly_fit_length(rop, op->length);
        _fq_poly_scalar_submul_fq(rop->coeffs, op->coeffs, op->length, x, ctx);
        _fq_poly_set_length(rop, FLINT_MAX(rop->length, op->length));
        _fq_poly_normalise(rop);
    }
}
