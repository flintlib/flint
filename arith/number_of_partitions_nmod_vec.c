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
arith_number_of_partitions_nmod_vec(mp_ptr res, slong len, nmod_t mod)
{
    mp_ptr tmp;
    mp_limb_t r;
    slong k, n;

    r = mod.n - UWORD(1);

    if (len < 1)
        return;

    tmp = _nmod_vec_init(len);
    _nmod_vec_zero(tmp, len);

    tmp[0] = UWORD(1);

    for (n = k = 1; n + 4*k + 2 < len; k += 2)
    {
        tmp[n] = r;
        tmp[n + k] = r;
        tmp[n + 3*k + 1] = UWORD(1);
        tmp[n + 4*k + 2] = UWORD(1);
        n += 6*k + 5;
    }

    if (n < len) tmp[n] = r;
    if (n + k < len) tmp[n + k] = r;
    if (n + 3*k + 1 < len) tmp[n + 3*k + 1] = WORD(1);

    _nmod_poly_inv_series(res, tmp, len, len, mod);

    _nmod_vec_clear(tmp);
}
