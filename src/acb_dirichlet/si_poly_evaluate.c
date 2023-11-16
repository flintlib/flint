/*
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* bsgs evaluation */
void
acb_dirichlet_si_poly_evaluate(acb_t res, slong * v, slong len, const acb_t z, slong prec)
{
    slong k, r, m;
    acb_t sq;
    acb_ptr zk;

    if (len < 3)
    {
        if (len == 0)
        {
            acb_zero(res);
        }
        else if (len == 1)
        {
            acb_set_si(res, v[0]);
        }
        else if (len == 2)
        {
            acb_mul_si(res, z, v[1], prec);
            acb_add_si(res, res, v[0], prec);
        }
        return;
    }

    m = n_sqrt(len) + 1;

    zk = _acb_vec_init(m + 1);
    _acb_vec_set_powers(zk, z, m + 1, prec);

    acb_init(sq);
    acb_zero(res);

    k = len - 1;
    r = k % m;
    for (; k >= 0; r = m - 1)
    {
        acb_dot_si(sq, NULL, 0, zk + r, -1, v + k, -1, r + 1, prec);
        k -= (r + 1);
        acb_mul(res, res, zk + m, prec);
        acb_add(res, res, sq, prec);
    }

    _acb_vec_clear(zk, m + 1);
    acb_clear(sq);
}
