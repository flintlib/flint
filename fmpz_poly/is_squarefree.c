/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

int _fmpz_poly_is_squarefree(const fmpz * poly, slong len)
{
    if (len < 3)
    {
        return 1;
    }
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
        int ans;
        fmpz * w = _fmpz_vec_init(2 * len);
        
        _fmpz_poly_derivative(w, poly, len);
        _fmpz_poly_gcd(w + len, poly, len, w, len - WORD(1));
        ans = _fmpz_vec_is_zero(w + len + 1, len - 2);
 
        _fmpz_vec_clear(w, 2 * len);
        return ans;
    }
}

int fmpz_poly_is_squarefree(const fmpz_poly_t poly)
{
    return _fmpz_poly_is_squarefree(poly->coeffs, poly->length);
}

