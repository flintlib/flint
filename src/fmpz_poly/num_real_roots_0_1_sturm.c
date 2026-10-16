/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

slong fmpz_poly_num_real_roots_0_1_sturm(const fmpz_poly_t pol)
{
    if (fmpz_poly_is_zero(pol))
        flint_throw(FLINT_ERROR, "(fmpz_poly_num_real_roots_0_1_sturm): zero polynomial\n");

    return _fmpz_poly_num_real_roots_sturm_bounded(pol->coeffs, pol->length, 1, WORD_MAX);
}
