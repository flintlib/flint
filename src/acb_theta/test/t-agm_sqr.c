/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_sqr....");
    fflush(stdout);
    
    flint_randinit(state);
    
    /* Test: agrees with agm_mul */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 500);
        slong bits = n_randint(state, 5);
        slong k;

        acb_ptr a, b, c;
        a = _acb_vec_init(n);
        b = _acb_vec_init(n);
        c = _acb_vec_init(n);

        for (k = 0; k < n; k++)
        {
            acb_randtest_precise(&a[k], state, prec, bits);
        }
        acb_theta_agm_sqr(b, a, g, prec);
        acb_theta_agm_mul(c, a, a, g, prec);

        if (!_acb_vec_overlaps(b, c, n))
        {
            flint_printf("FAIL\n\n");
            flint_abort();
        }

        _acb_vec_clear(a, n);
        _acb_vec_clear(b, n);
        _acb_vec_clear(c, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
