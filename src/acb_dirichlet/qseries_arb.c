/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* sum a_k x^(k^2), assume a[0] = 0 */
void
acb_dirichlet_qseries_arb(acb_t res, acb_srcptr a, const arb_t x, slong len, slong prec)
{
    slong k;
    arb_t xk2, dx, x2;

    arb_init(xk2);
    arb_init(dx);
    arb_init(x2);

    arb_set(dx, x);
    arb_set(xk2, dx);
    arb_mul(x2, dx, dx, prec);

    acb_mul_arb(res, a + 1, xk2, prec);

    /* TODO: reduce prec */
    for (k = 2; k < len; k++)
    {
        arb_mul(dx, dx, x2, prec);
        arb_mul(xk2, xk2, dx, prec);
        acb_addmul_arb(res, a + k, xk2, prec);
    }

    arb_clear(xk2);
    arb_clear(x2);
    arb_clear(dx);
}
