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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

int _fmpq_poly_is_squarefree(const fmpz * poly, const fmpz_t den, len_t len)
{
    if (len < 3)
        return 1;
    else if (len == 3)
    {
        int ans;
        fmpz_t lhs, rhs;
        fmpz_init(lhs);
        fmpz_init(rhs);
        
        fmpz_mul(lhs, poly + 1, poly + 1);
        fmpz_mul(rhs, poly, poly + 2);
        fmpz_mul_ui(rhs, rhs, 4);

        ans = !fmpz_equal(lhs, rhs);
        fmpz_clear(lhs);
        fmpz_clear(rhs);
        return ans;
    }
    else
    {
        len_t gdeg;
        fmpz * w = _fmpz_vec_init(2 * len);
        
        _fmpz_poly_derivative(w, poly, len);
        _fmpz_poly_gcd(w + len, poly, len, w, len - 1L);
        
        for (gdeg = len - 2L; w[gdeg] == 0L; gdeg--) ;
        
        _fmpz_vec_clear(w, 2 * len);
        return (gdeg == 0);
    }
}

int fmpq_poly_is_squarefree(const fmpq_poly_t poly)
{
    return _fmpq_poly_is_squarefree(poly->coeffs, poly->den, poly->length);
}

