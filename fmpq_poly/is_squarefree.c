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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

int _fmpq_poly_is_squarefree(const fmpz * poly, const fmpz_t den, long len)
{
    if (len < 3)
        return 1;
    else if (len == 3)
    {
        int ans;
        fmpz_t lhs, rhs;
        fmpz_init(lhs);
        fmpz_init(rhs);
        fmpz_mul(lhs, poly + 2, poly);
        fmpz_mul_ui(lhs, lhs, 4UL);
        fmpz_mul(rhs, poly + 1, poly + 1);
        ans = !fmpz_equal(lhs, rhs);
        fmpz_clear(lhs);
        fmpz_clear(rhs);
        return ans;
    }
    else
    {
        long gdeg;
        fmpz * der, * gcd;
        fmpz * Z = _fmpz_vec_init(2 * (len - 1));
        der = Z;
        gcd = Z + len - 1;
        _fmpz_poly_derivative(der, poly, len);
        _fmpz_poly_gcd(gcd, poly, len, der, len - 1);
        for (gdeg = len - 2; fmpz_is_zero(gcd + gdeg); gdeg--) ;
        _fmpz_vec_clear(Z, 2 * (len - 1));
        return gdeg < 1;
    }
}

int fmpq_poly_is_squarefree(const fmpq_poly_t poly)
{
    return _fmpq_poly_is_squarefree(poly->coeffs, poly->den, poly->length);
}

