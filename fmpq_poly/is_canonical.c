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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

int _fmpq_poly_is_canonical(const fmpz * poly, const fmpz_t den, len_t len)
{
    if (len)
    {
        int ans;
        fmpz_t c;

        if (fmpz_is_zero(poly + len - 1))
            return 0;

        if (fmpz_sgn(den) < 0)
            return 0;

        fmpz_init(c);
        _fmpz_poly_content(c, poly, len);
        fmpz_gcd(c, c, den);
        ans = (*c == 1L);
        fmpz_clear(c);

        return ans;
    }
    else
    {
        return (*den == 1L);
    }
}

int fmpq_poly_is_canonical(const fmpq_poly_t poly)
{
    return _fmpq_poly_is_canonical(poly->coeffs, poly->den, poly->length);
}

