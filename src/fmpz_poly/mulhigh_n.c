/*
    Copyright (C) 2008, 2009 William Hart

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

void
fmpz_poly_mulhigh_n(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
{
    mp_size_t limbs1 = _fmpz_vec_max_limbs(poly1->coeffs, poly1->length);
    mp_size_t limbs2 = _fmpz_vec_max_limbs(poly2->coeffs, poly2->length);
    mp_size_t len1 = poly1->length;
    mp_size_t len2 = poly2->length;
    mp_size_t limbsx = FLINT_MAX(limbs1, limbs2);

    if (n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (n < 4)
    {
        fmpz_poly_mulhigh_classical(res, poly1, poly2, n - 1);
        return;
    }

    if ((limbsx > 4) && (n < 16) && poly1->length <= n && poly2->length <= n)
        fmpz_poly_mulhigh_karatsuba_n(res, poly1, poly2, n);
    else if (limbs1 + limbs2 <= 8)
        fmpz_poly_mul_KS(res, poly1, poly2);
    else if ((limbs1+limbs2)/2048 > len1 + len2)
        fmpz_poly_mul_KS(res, poly1, poly2);
    else if ((limbs1 + limbs2)*FLINT_BITS*4 < len1 + len2)
       fmpz_poly_mul_KS(res, poly1, poly2);
    else
       fmpz_poly_mul_SS(res, poly1, poly2);
}
