/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

static void
bsplit_nmod(mp_ptr R, mp_ptr Q, slong a, slong b, nmod_t mod)
{
    if (b - a == 1)
    {
        R[0] = 0;
        R[1] = 1;
        Q[0] = 1;
        Q[1] = nmod_neg(nmod_set_ui(a, mod), mod);
    }
    else
    {
        slong m, len1, len2;
        mp_ptr R1, R2, Q1, Q2;
        m = a + (b - a) / 2;

        len1 = (m - a) + 1;
        len2 = (b - m) + 1;

        R1 = _nmod_vec_init(2 * len1 + 2 * len2);
        Q1 = R1 + len1;
        R2 = Q1 + len1;
        Q2 = R2 + len2;

        bsplit_nmod(R1, Q1, a, m, mod);
        bsplit_nmod(R2, Q2, m, b, mod);

        _nmod_poly_mul(R, Q2, len2, R1, len1, mod);
        _nmod_vec_add(R + (m - a), R + (m - a), R2, len2, mod);
        _nmod_poly_mul(Q, Q2, len2, Q1, len1, mod);

        _nmod_vec_clear(R1);
    }
}

void
arith_bell_number_nmod_vec_ogf(mp_ptr res, slong len, nmod_t mod)
{
    mp_ptr R, Q;

    if (len <= 2 || mod.n == 1)
    {
        arith_bell_number_nmod_vec_recursive(res, len, mod);
        return;
    }

    R = _nmod_vec_init(len + 1);
    Q = _nmod_vec_init(len + 1);

    bsplit_nmod(R, Q, 1, len, mod);
    _nmod_poly_div_series(res, R, len, Q, len, len, mod);
    res[0] = 1;

    _nmod_vec_clear(R);
    _nmod_vec_clear(Q);
}
