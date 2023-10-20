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

    flint_printf("precomp_set....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: can clear and init again; values are all 1 if input is all zero */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong nb = n_randint(state, 3);
        slong prec = ACB_THETA_LOW_PREC;
        slong mag_bits = n_randint(state, 2);
        acb_theta_eld_t E;
        acb_theta_precomp_t D;
        acb_mat_t tau;
        arb_mat_t C;
        acb_ptr zs;
        arb_t x;
        arf_t R2;
        arb_ptr v;
        slong k, j;
        int res = 1;

        acb_theta_eld_init(E, g, g);
        acb_theta_precomp_init(D, nb, g);
        acb_mat_init(tau, g, g);
        arb_mat_init(C, g, g);
        arf_init(R2);
        v = _arb_vec_init(g);
        zs = _acb_vec_init(g * nb);
        arb_init(x);

        acb_theta_precomp_clear(D);
        acb_theta_precompo_init(D, nb, g);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_theta_eld_cho(C, tau, prec);
        arb_randtest_positive(x, state, prec, mag_bits);
        arf_set(R2, arb_midref(x));
        arf_mul_si(R2, R2, 1 + n_randint(state, 10), prec, ARF_RND_UP);
        for (k = 0; k < g; k++)
        {
            arb_randtest_precise(&v[k], state, prec, mag_bits);
        }
        acb_theta_eld_fill(E, C, R2, v);
        acb_mat_zero(tau);

        acb_theta_precomp_set(D, zs, tau, E, prec);

        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                if (!acb_is_one(acb_mat_entry(acb_theta_precomp_exp_mat(D), j, k)))
                {
                    res = 0;
                }
            }
        }

        for (k = 0; k < g; k++)
        {
            for (j = 0; j <= acb_theta_eld_box(E, k); j++)
            {
                if (!acb_is_one(acb_theta_precomp_sqr_pow(D, k, j)))
                {
                    res = 0;
                }
            }
        }

        for (k = 0; k < nb; k++)
        {
            for (j = 0; j < g; j++)
            {
                if (!acb_is_one(acb_theta_precomp_exp_z(D, k, j)))
                {
                    res = 0;
                }
            }
        }

        if (!res)
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_theta_eld_clear(E);
        acb_theta_precomp_clear(D);
        acb_mat_clear(tau);
        arb_mat_clear(C);
        arf_clear(R2);
        _arb_vec_clear(v, g);
        _acb_vec_clear(zs, nb * g);
        arb_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
