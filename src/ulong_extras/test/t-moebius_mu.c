/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

#ifndef check
#define check check
void check(mp_limb_t n, int s1, int s2)
{
    if (s1 != s2)
    {
        flint_printf("FAIL:\n");
        flint_printf("%wu: got %d instead of %d\n", n, s1, s2);
        fflush(stdout);
        flint_abort();
    }
}
#endif

TEST_FUNCTION_START(n_moebius_mu, state)
{
    int n, k, s;
    int * mu;

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
        fflush(stdout);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
