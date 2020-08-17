/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_rank(const fmpz_mat_t A)
{
    fmpz_mat_t tmp;
    fmpz_t den;
    slong rank;

    if (fmpz_mat_is_empty(A))
        return 0;

    fmpz_mat_init_set(tmp, A);
    fmpz_init(den);

    if (FLINT_ABS(fmpz_mat_max_bits(tmp)) <= (FLINT_BITS - 4)/2)
        rank = fmpz_mat_rank_small_inplace(tmp);
    else if (FLINT_MIN(tmp->r, tmp->c) < 25)
        rank = fmpz_mat_fflu(tmp, den, NULL, tmp, 0);
    else
        rank = fmpz_mat_rref(tmp, den, tmp);

    fmpz_mat_clear(tmp);
    fmpz_clear(den);

    return rank;
}
