/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void
arith_number_of_partitions_vec(fmpz * res, slong len)
{
    fmpz * tmp;
    slong k, n;

    if (len < 1)
        return;

    tmp = _fmpz_vec_init(len);

    tmp[0] = WORD(1);

    for (n = k = 1; n + 4*k + 2 < len; k += 2)
    {
        tmp[n] = WORD(-1);
        tmp[n + k] = WORD(-1);
        tmp[n + 3*k + 1] = WORD(1);
        tmp[n + 4*k + 2] = WORD(1);
        n += 6*k + 5;
    }

    if (n < len) tmp[n] = WORD(-1);
    if (n + k < len) tmp[n + k] = WORD(-1);
    if (n + 3*k + 1 < len) tmp[n + 3*k + 1] = WORD(1);

    _fmpz_poly_inv_series(res, tmp, len, len);
     
    _fmpz_vec_clear(tmp, len);
}
