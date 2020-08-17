/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"


void check(mp_limb_t n, int mu1, int mu2)
{
    if (mu1 != mu2)
    {
        flint_printf("FAIL:\n");
        flint_printf("mu(%wu): %d != %d\n", n, mu1, mu2); 
        abort();
    }
}

int main(void)
{
    int n, k, s;
    int * mu;

    FLINT_TEST_INIT(state);
    
    flint_printf("moebius_mu....");
    fflush(stdout);

    check(0, n_moebius_mu(0), 0);
    check(1, n_moebius_mu(1), 1);

    for (n = 1; n < 100; n++)
    {
        mu = flint_malloc(sizeof(int) * n);
        n_moebius_mu_vec(mu, n);
        for (k = 0; k < n; k++)
            check(k, mu[k], n_moebius_mu(k));
        flint_free(mu);
    }

    mu = flint_malloc(sizeof(int) * 10000);
    n_moebius_mu_vec(mu, 10000);
    for (k = 0; k < n; k++)
        check(k, mu[k], n_moebius_mu(k));
    flint_free(mu);

    check(10000, n_moebius_mu(10000), 0);
    check(10001, n_moebius_mu(10001), 1);
    check(10002, n_moebius_mu(10002), -1);
    check(10003, n_moebius_mu(10003), 1);
    check(10004, n_moebius_mu(10004), 0);
    check(10005, n_moebius_mu(10005), 1);
    check(10006, n_moebius_mu(10006), 1);
    check(10007, n_moebius_mu(10007), -1);
    check(10008, n_moebius_mu(10008), 0);
    check(10009, n_moebius_mu(10009), -1);
    check(10010, n_moebius_mu(10010), -1);

    s = 0;
    for (k = 0; k <= 10000; k++)
        s += n_moebius_mu(k);

    if (s != -23)
    {
        flint_printf("FAIL:\n");
        flint_printf("expected mu(k), k <= 10000 to sum to %d (got %d)\n", -23, s);
        abort();
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
