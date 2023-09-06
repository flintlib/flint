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

    flint_printf("jet_fourier....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: doing it twice gives rescaled & reordered vector */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong prec = ACB_THETA_LOW_PREC + n_randint(state, 1000);
        slong ord = n_randint(state, 4);
        slong g = 1 + n_randint(state, 4);
        slong b = ord + 1;
        slong nb = n_pow(b, g);
        acb_ptr val, tf, test;
        acb_t c;
        slong k, kk, j, i;

        val = _acb_vec_init(nb);
        tf = _acb_vec_init(nb);
        test = _acb_vec_init(nb);
        acb_init(c);

        for (k = 0; k < nb; k++)
        {
            acb_urandom(&val[k], state, prec);
        }

        acb_theta_jet_fourier(tf, val, ord, g, prec);
        acb_theta_jet_fourier(tf, tf, ord, g, prec);

        /* Set test to rescaled & reordered vector */
        for (k = 0; k < nb; k++)
        {
            kk = k;
            j = 0;
            for (i = 0; i < g; i++)
            {
                j += ((b - (kk % b)) % b) * n_pow(b, i);
                kk = kk / b;
            }
            acb_set(&test[j], &val[k]);
        }
        acb_set_si(c, b);
        acb_pow_ui(c, c, g, prec);
        _acb_vec_scalar_mul(test, test, nb, c, prec);

        if (!_acb_vec_overlaps(test, tf, nb))
        {
            flint_printf("FAIL (double transform)\n");
            flint_printf("g = %wd, ord = %wd, values:\n", g, ord);
            _acb_vec_printd(val, nb, 5);
            flint_printf("transform:\n");
            _acb_vec_printd(tf, nb, 5);
            flint_printf("rescaled:\n");
            _acb_vec_printd(test, nb, 5);
            flint_abort();
        }

        _acb_vec_clear(val, nb);
        _acb_vec_clear(tf, nb);
        _acb_vec_clear(test, nb);
        acb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
