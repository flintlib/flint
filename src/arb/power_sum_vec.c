/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "bernoulli.h"
#include "arb_poly.h"

/* todo: don't use exact bernoulli numbers for large len */
/* todo: output exact integers when precise enough */
void
arb_power_sum_vec(arb_ptr res, const arb_t a, const arb_t b, slong len, slong prec)
{
    arb_ptr t, u, v;
    slong k;

    if (len < 1)
        return;

    t = _arb_vec_init(len + 1);
    u = _arb_vec_init(len + 1);
    v = _arb_vec_init(len + 1);

    /* exp(ax), exp(bx) */
    arb_set(t + 1, a);
    arb_set(u + 1, b);
    _arb_poly_exp_series(t, t, 2, len + 1, prec);
    _arb_poly_exp_series(u, u, 2, len + 1, prec);
    _arb_vec_sub(t, u, t, len + 1, prec);

    /* x/(exp(x)-1) */
    BERNOULLI_ENSURE_CACHED(len + 1);
    for (k = 0; k <= len; k++)
        arb_set_fmpq(u + k, bernoulli_cache + k, prec);
    _arb_poly_borel_transform(u, u, len + 1, prec);

    _arb_poly_mullow(v, t, len + 1, u, len + 1, len + 1, prec);
    _arb_poly_inv_borel_transform(v, v, len + 1, prec);

    for (k = 0; k < len; k++)
        arb_div_ui(res + k, v + k + 1, k + 1, prec);

    _arb_vec_clear(t, len + 1);
    _arb_vec_clear(u, len + 1);
    _arb_vec_clear(v, len + 1);
}

