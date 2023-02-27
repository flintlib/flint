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
_fmpz_poly_mulhigh(fmpz * res, const fmpz * poly1, slong len1,
                                   const fmpz * poly2, slong len2, slong start)
{
    mp_size_t limbs1 = _fmpz_vec_max_limbs(poly1, len1);
    mp_size_t limbs2 = _fmpz_vec_max_limbs(poly2, len2);
    mp_size_t limbsx = FLINT_MAX(limbs1, limbs2);

    if (start < 5)
    {
        _fmpz_poly_mulhigh_classical(res, poly1, len1, poly2, len2, start);
        return;
    }

    if (limbsx > 4 && start < 17 && len1 == start + 1 && len2 == start + 1)
        _fmpz_poly_mulhigh_karatsuba_n(res, poly1, poly2, start + 1);
    else if (limbs1 + limbs2 <= 8)
        _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
    else if ((limbs1 + limbs2)/2048 > len1 + len2)
        _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
    else if ((limbs1 + limbs2)*FLINT_BITS*4 < len1 + len2)
       _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
    else
       _fmpz_poly_mul_SS(res, poly1, len1, poly2, len2);
}
