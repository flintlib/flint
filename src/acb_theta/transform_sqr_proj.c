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
acb_theta_transform_sqr_proj(acb_ptr res, acb_srcptr th2, const fmpz_mat_t mat, slong prec)
{
    acb_ptr aux;
    slong g = sp2gz_dim(mat);
    ulong n2 = 1 << (2 * g);
    ulong ab;
    ulong image_ab;
    fmpz_t eps;
    acb_t c;

    aux = _acb_vec_init(n2);
    fmpz_init(eps);
    acb_init(c);

    for (ab = 0; ab < n2; ab++)
    {
        image_ab = acb_theta_transform_char(eps, ab, mat);
        acb_unit_root(c, 4, prec); /* 8 for theta values, 4 for squares */
        acb_pow_fmpz(c, c, eps, prec);
        acb_mul(c, c, &th2[image_ab], prec);
        acb_set(&aux[ab], c);
    }

    _acb_vec_set(res, aux, n2);

    _acb_vec_clear(aux, n2);
    fmpz_clear(eps);
    acb_clear(c);
}
