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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"


void
_fmpz_poly_div_root(fmpz * Q, const fmpz * A, long len, const fmpz_t c)
{
    fmpz_t r, t;
    long i;

    if (len < 2)
        return;

    fmpz_init(r);
    fmpz_init(t);
    fmpz_set(r, A + len - 1);

    for (i = len - 2; i > 0; i--)
    {
        fmpz_mul(t, r, c);
        fmpz_add(t, t, A + i);
        fmpz_swap(Q + i, r);
        fmpz_swap(r, t);
    }

    fmpz_swap(Q, r);
    fmpz_clear(r);
    fmpz_clear(t);
}

void
fmpz_poly_div_root(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_t c)
{
    long len = A->length;

    if (len <= 1)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (fmpz_is_zero(c))
    {
        fmpz_poly_shift_right(Q, A, 1);
        return;
    }

    fmpz_poly_fit_length(Q, len - 1);
    _fmpz_poly_div_root(Q->coeffs, A->coeffs, len, c);
    _fmpz_poly_set_length(Q, len - 1);
}
