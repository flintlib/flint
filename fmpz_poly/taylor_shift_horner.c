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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_taylor_shift_horner(fmpz * poly, const fmpz_t c, len_t n)
{
    len_t i, j;

    if (*c == 1L)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                fmpz_add(poly + j, poly + j, poly + j + 1);
    }
    else if (*c == -1L)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                fmpz_sub(poly + j, poly + j, poly + j + 1);
    }
    else if (*c != 0L)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                fmpz_addmul(poly + j, poly + j + 1, c);
    }
}

void
fmpz_poly_taylor_shift_horner(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_taylor_shift_horner(g->coeffs, c, g->length);
}
