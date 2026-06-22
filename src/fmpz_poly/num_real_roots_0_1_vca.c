/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

slong fmpz_poly_num_real_roots_0_1_vca(const fmpz_poly_t pol)
{
    slong n_exact = 0;
    slong n_interval = 0;

    _fmpz_poly_isolate_real_roots_0_1_vca(NULL, &n_exact, NULL, NULL, &n_interval, pol->coeffs, pol->length);

    return n_exact + n_interval;
}
