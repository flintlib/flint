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
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

void
fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, const nmod_poly_t poly)
{
    long len = poly->length;

    if (len == 0)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        long i;
        fmpz_poly_fit_length(res, len);
        for (i = 0; i < len; i++)
            fmpz_set_ui(res->coeffs + i, poly->coeffs[i]);
        _fmpz_poly_set_length(res, len);
    }
}
