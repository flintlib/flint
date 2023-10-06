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

    flint_printf("eld_border....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: border points are not contained in the ellipsoid */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = ACB_THETA_LOW_PREC;
        slong mag_bits = n_randint(state, 2);
        acb_theta_eld_t E;
        arb_mat_t C;
        arb_t x;
        arf_t R2;
        arb_ptr v;
        slong k;
        slong *all_pts;

        acb_theta_eld_init(E, g, g);
        arb_mat_init(C, g, g);
        arf_init(R2);
        v = _arb_vec_init(g);
        arb_init(x);

        arb_mat_randtest_cho(C, state, prec, mag_bits);
        arb_mat_transpose(C, C);
        arb_randtest_positive(x, state, prec, mag_bits);
        arf_set(R2, arb_midref(x));
        arf_mul_si(R2, R2, 1 + n_randint(state, 10), prec, ARF_RND_UP);
        for (k = 0; k < g; k++)
        {
            arb_randtest_precise(&v[k], state, prec, mag_bits);
        }

        acb_theta_eld_fill(E, C, R2, v);
        all_pts = flint_malloc(acb_theta_eld_nb_border(E) * g * sizeof(slong));
        acb_theta_eld_border(all_pts, E);

        for (k = 0; k < acb_theta_eld_nb_border(E); k++)
        {
            for (j = 0; j < g; j++)
            {
                if (acb_theta_eld_contains(E, all_pts + k * g))
                {
                    flint_printf("FAIL: point inside ellipsoid\n");
                    flint_printf("\n");
                    flint_abort();
                }
            }
        }

        acb_theta_eld_clear(E);
        arb_mat_clear(C);
        arf_clear(R2);
        _arb_vec_clear(v, g);
        flint_free(all_pts);
        arb_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
