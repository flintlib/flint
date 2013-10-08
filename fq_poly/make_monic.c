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

    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include "fq_poly.h"

void
_fq_poly_make_monic(fq_struct *rop, 
                    const fq_struct *op, long length, const fq_ctx_t ctx)
{
    fq_t inv;
    fq_init(inv);
    fq_inv(inv, &op[length - 1], ctx);
    _fq_poly_scalar_mul_fq(rop, op, length, inv, ctx);
    fq_clear(inv);
}

void
fq_poly_make_monic(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
{
    if (op->length == 0)
    {
        fq_poly_zero(rop);
    }

    fq_poly_fit_length(rop, op->length);
    _fq_poly_make_monic(rop->coeffs, 
                            op->coeffs, op->length, ctx);
    _fq_poly_set_length(rop, op->length);
}

