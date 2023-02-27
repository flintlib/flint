/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"


slong
fmpz_poly_mat_rank(const fmpz_poly_mat_t A)
{
    fmpz_poly_mat_t tmp;
    fmpz_poly_t den;
    slong rank;

    if (fmpz_poly_mat_is_empty(A))
        return 0;

    fmpz_poly_mat_init_set(tmp, A);
    fmpz_poly_init(den);
    rank = fmpz_poly_mat_fflu(tmp, den, NULL, tmp, 0);
    fmpz_poly_mat_clear(tmp);
    fmpz_poly_clear(den);
    return rank;
}
