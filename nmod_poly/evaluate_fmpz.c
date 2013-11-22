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

    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

/* TODO: Add a divconquer method */

void
_nmod_poly_evaluate_fmpz(fmpz_t rop, const mp_srcptr poly, const slong len, const fmpz_t c)
{
    fmpz_t t;
    slong m;

    if (len == 0)
    {
        fmpz_zero(rop);
        return;
    }

    if (len == 1 || c == 0)
    {
        fmpz_set_ui(rop, poly[0]);
        return;
    }

    m = len - 1;
    

    fmpz_init(t);
    fmpz_set_ui(rop, poly[m]);
    m--;

    for ( ; m >= 0; m--)
    {
        fmpz_mul(t, rop, c);
        fmpz_add_ui(rop, t, poly[m]);
    }
    fmpz_clear(t);
}

void
nmod_poly_evaluate_fmpz(fmpz_t rop, const nmod_poly_t poly, const fmpz_t c)
{
    _nmod_poly_evaluate_fmpz(rop, poly->coeffs, poly->length, c);
}

