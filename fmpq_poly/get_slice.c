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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_get_slice(fmpq_poly_t rop, const fmpq_poly_t op, long i, long j)
{
    i = FLINT_MAX(i, 0);
    j = FLINT_MIN(j, op->length);

    if (i < j)
    {
        long k;

        if (rop == op)
        {
            for (k = 0; k < i; k++)
                fmpz_zero(rop->coeffs + k);
            for (k = j; k < rop->length; k++)
                fmpz_zero(rop->coeffs + k);
            fmpq_poly_canonicalise(rop);
        }
        else
        {
            fmpq_poly_fit_length(rop, j);
            _fmpq_poly_set_length(rop, j);

            _fmpz_vec_set(rop->coeffs + i, op->coeffs + i, j - i);
            fmpz_set(rop->den, op->den);
            fmpq_poly_canonicalise(rop);
        }
    }
    else
    {
        fmpq_poly_zero(rop);
    }
}

