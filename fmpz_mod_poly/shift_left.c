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

    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009, 2008 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly, len_t len, len_t n)
{
    len_t i;

    /* Copy in reverse to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = len; i--; )
            fmpz_set(res + n + i, poly + i);
    }
    else
    {
        for (i = len; i--; )
            fmpz_swap(res + n + i, res + i);
    }

    for (i = 0; i < n; i++)
        fmpz_zero(res + i);
}

void
fmpz_mod_poly_shift_left(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, len_t n)
{
    if (n == 0)
    {
        fmpz_mod_poly_set(res, poly);
        return;
    }

    if (poly->length == 0)
    {
        fmpz_mod_poly_zero(res);
        return;
    }

    fmpz_mod_poly_fit_length(res, poly->length + n);
    _fmpz_mod_poly_shift_left(res->coeffs, poly->coeffs, poly->length, n);
    _fmpz_mod_poly_set_length(res, poly->length + n);
}

