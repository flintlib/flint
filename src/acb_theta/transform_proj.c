/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_transform_proj(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th, int sqr, slong prec)
{
    slong g = sp2gz_dim(mat);
    ulong n2 = 1 << (2 * g);
    slong k = (sqr ? 4 : 8);
    acb_ptr aux;
    ulong ab;
    ulong image_ab;
    fmpz_t eps;
    acb_t c;

    aux = _acb_vec_init(n2);
    fmpz_init(eps);
    acb_init(c);

    for (ab = 0; ab < n2; ab++)
    {
        image_ab = acb_theta_transform_char(eps, mat, ab);
        acb_unit_root(c, k, prec);
        acb_pow_fmpz(c, c, eps, prec);
        acb_mul(c, c, &th[image_ab], prec);
        acb_set(&aux[ab], c);
    }

    _acb_vec_set(res, aux, n2);

    _acb_vec_clear(aux, n2);
    fmpz_clear(eps);
    acb_clear(c);
}
