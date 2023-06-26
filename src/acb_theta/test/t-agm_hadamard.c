/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_hadamard....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: twice Hadamard should be multiplication by 2^g */
    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 5);
        slong prec = 20 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 4);
        acb_ptr s;
        acb_ptr r;
        acb_ptr test;
        slong n = 1 << g;
        slong k;
        int res;

        s = _acb_vec_init(n);
        r = _acb_vec_init(n);
        test = _acb_vec_init(n);

        for (k = 0; k < n; k++)
            acb_randtest_precise(&s[k], state, prec, mag_bits);
        acb_theta_agm_hadamard(r, s, g, prec);
        acb_theta_agm_hadamard(test, r, g, prec);
        _acb_vec_scalar_mul_2exp_si(test, test, n, -g);

        res = 1;
        for (k = 0; k < n; k++)
        {
            if (!acb_contains(&test[k], &s[k]))
                res = 0;
        }
        if (!res)
        {
            flint_printf("FAIL (overlap):\n");
            for (k = 0; k < n; k++)
            {
                acb_printd(&s[k], 10);
                flint_printf("\n");
                acb_printd(&test[k], 10);
                flint_printf("\n");
            }
            fflush(stdout);
            flint_abort();
        }

        _acb_vec_clear(s, n);
        _acb_vec_clear(r, n);
        _acb_vec_clear(test, n);
    }


    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
