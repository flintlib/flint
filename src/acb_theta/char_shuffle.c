/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_theta_char_shuffle(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th, int sqr, slong prec)
{
    slong g = sp2gz_dim(mat);
    ulong n2 = 1 << (2 * g);
    slong k = (sqr ? 4 : 8);
    acb_ptr aux, units;
    ulong ab;
    ulong * chars;
    slong * es;

    aux = _acb_vec_init(n2);
    units = _acb_vec_init(k);
    chars = flint_malloc(n2 * sizeof(ulong));
    es = flint_malloc(n2 * sizeof(slong));

    acb_theta_char_table(chars, es, mat, -1);
    _acb_vec_unit_roots(units, k, k, prec);

    for (ab = 0; ab < n2; ab++)
    {
        acb_mul(&aux[ab], &units[es[ab] % k], &th[chars[ab]], prec);
    }
    _acb_vec_set(res, aux, n2);

    _acb_vec_clear(aux, n2);
    _acb_vec_clear(units, k);
    flint_free(es);
    flint_free(chars);
}
