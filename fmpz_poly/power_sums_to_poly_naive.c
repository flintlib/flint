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

    Copyright (C) 2016 Vincent Delecroix

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_power_sums_to_poly_naive(fmpz * res, const fmpz * poly,
                                       slong len)
{
    slong i, k;
    slong d = fmpz_get_ui(poly);

    fmpz_set_ui(res + d, 1);
    for (k = 1; k < FLINT_MIN(d + 1, len); ++k)
    {
        fmpz_set(res + d - k, poly + k);
        for (i = 1; i < k; ++i)
            fmpz_addmul(res + d - k, res + d - k + i, poly + i);
        fmpz_divexact_si(res + d - k, res + d - k, k);
        fmpz_neg(res + d - k, res + d - k);
    }
    for (k = len; k <= d; ++k)
    {
        fmpz_set(res + d - k, poly + k);
        for (i = 1; i < len; ++i)
            fmpz_addmul(res + d - k, res + d - k + i, poly + i);
        fmpz_divexact_si(res + d - k, res + d - k, k);
        fmpz_neg(res + d - k, res + d - k);
    }

}
