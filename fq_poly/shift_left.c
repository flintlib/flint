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

void _fq_poly_shift_left(fq_struct *rop, const fq_struct *op, long len, long n)
{
    long i;

    if (rop != op)
    {
        for (i = len; i--; )
            fq_set(rop + n + i, op + i);
    }
    else
    {
        for (i = len; i--; )
            fq_swap(rop + n + i, rop + i);
    }

    for (i = 0; i < n; i++)
        fq_zero(rop + i);
}

void fq_poly_shift_left(fq_poly_t rop, const fq_poly_t op, long n)
{
    if (n == 0)
    {
        fq_poly_set(rop, op);
    }
    else if (fq_poly_is_zero(op))
    {
        fq_poly_zero(rop);
    }
    else
    {
        fq_poly_fit_length(rop, op->length + n);
        _fq_poly_shift_left(rop->coeffs, op->coeffs, op->length, n);
        _fq_poly_set_length(rop, op->length + n);
    }
}
