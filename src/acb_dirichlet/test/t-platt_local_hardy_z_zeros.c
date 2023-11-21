/*
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_platt_local_hardy_z_zeros, state)
{
    /* Check a specific combination of parameter values that is relatively fast
     * to evaluate and that has relatively tight bounds. */
    slong A, B, K, sigma_grid, Ns_max, sigma_interp;
    arb_t h, H;
    fmpz_t J, T, n;
    arb_ptr pa, pb;
    slong count, i;
    slong maxcount = 50;
    slong prec = 128;

    arb_init(h);
    arb_init(H);
    fmpz_init(J);
    fmpz_init(T);
    fmpz_init(n);
    pa = _arb_vec_init(maxcount);
    pb = _arb_vec_init(maxcount);

    fmpz_set_si(n, 10142);

    /* parameters related to the location/resolution/width of the grid */
    fmpz_set_si(T, 10000);
    A = 8;
    B = 128;

    /* tuning parameters for the evaluation of grid points */
    fmpz_set_si(J, 1000);
    K = 30;
    sigma_grid = 63;
    arb_set_d(h, 4.5);

    /* tuning parameters for interpolation on the grid */
    Ns_max = 200;
    sigma_interp = 21;
    arb_one(H);

    count = _acb_dirichlet_platt_local_hardy_z_zeros(pa, n, maxcount,
            T, A, B, h, J, K, sigma_grid, Ns_max, H, sigma_interp, prec);
    acb_dirichlet_hardy_z_zeros(pb, n, count, prec);
    if (count != maxcount)
    {
        flint_printf("FAIL: not enough zeros were isolated\n\n");
        flint_printf("count = %wd  maxcount = %wd\n\n", count, maxcount);
        flint_abort();
    }

    for (i = 0; i < count; i++)
    {
        if (!arb_overlaps(pa+i, pb+i))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("observed[%wd] = ", i);
            arb_printd(pa+i, 20); flint_printf("\n\n");
            flint_printf("expected[%wd] = ", i);
            arb_printd(pb+i, 20); flint_printf("\n\n");
            flint_abort();
        }
    }

    arb_clear(h);
    arb_clear(H);
    fmpz_clear(J);
    fmpz_clear(T);
    fmpz_clear(n);
    _arb_vec_clear(pa, maxcount);
    _arb_vec_clear(pb, maxcount);

    TEST_FUNCTION_END(state);
}
