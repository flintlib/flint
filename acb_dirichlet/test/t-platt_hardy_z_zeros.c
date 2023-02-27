/*
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

int main()
{
    fmpz_t n;
    arb_ptr pa, pb;
    slong count, i;
    slong maxcount = 50;
    slong prec = 64;

    flint_printf("platt_hardy_z_zeros....");
    fflush(stdout);

    fmpz_init(n);
    pa = _arb_vec_init(maxcount);
    pb = _arb_vec_init(maxcount);

    fmpz_set_si(n, 10000);

    count = acb_dirichlet_platt_hardy_z_zeros(pa, n, maxcount, prec);
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

    fmpz_clear(n);
    _arb_vec_clear(pa, maxcount);
    _arb_vec_clear(pb, maxcount);

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
