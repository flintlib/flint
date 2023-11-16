/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_modular.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

#define EPS 1e-13
#define NUM_DERIVS 4
#define NUM_TESTS 7

const double k_testdata[NUM_TESTS][10] = {
    {0.0, 0.0, 1.5707963267948966192, 0.0, 0.39269908169872415481, 0.0,
        0.22089323345553233708, 0.0, 0.15339807878856412297, 0.0},
    {0.5, 0.0, 1.8540746773013719184, 0.0, 0.84721308479397908661, 0.0,
        0.92703733865068595922, 0.0, 1.2708196271909686299, 0.0},
    {-1.0, 0.0, 1.3110287771460599052, 0.0, 0.17798966494456595038, 0.0,
        0.051552950136795718707, 0.0, 0.018179887959689603011, 0.0},
    {2.0, 0.0, 1.3110287771460599052, -1.3110287771460599052,
        -0.47752472362846400224, 0.17798966494456595038, 0.2762042441497192576,
        -0.051552950136795718707, -0.18666835846938225718, 0.018179887959689603011},
    {-3.0, 0.0, 1.0782578237498216177, 0.0, 0.078788301660931975886, 0.0,
        0.011748068987044517782, 0.0, 0.0021065590680576326689, 0.0},
    {1.0, 1.0, 1.5092369540512728293, 0.62514641520269688427,
        -0.079689518666051625811, 0.40903679001382547524, -0.23159955416582020342,
        -0.028627375924621981072, 0.014191044435751759097, -0.16030448214657194629},
    {-2.0, -3.0, 1.0408718798817036, -0.24497111630480680352, 0.04401149835588265436,
        -0.060184042863324675054, -0.0012755513109907959184,
        -0.01044301570409968822, -0.0013811810360989366762, -0.0011248246747562196271}
};

TEST_FUNCTION_START(acb_modular_elliptic_k, state)
{
    slong iter;

    /* check particular values against table */
    {
        acb_t z, t;
        acb_ptr w1;
        slong i, j, prec, cnj;

        acb_init(z);
        acb_init(t);
        w1 = _acb_vec_init(NUM_DERIVS);

        for (prec = 32; prec <= 512; prec *= 4)
        {
            for (i = 0; i < NUM_TESTS; i++)
            {
                for (cnj = 0; cnj < 2; cnj++)
                {
                    if (cnj == 1 && k_testdata[i][0] > 1 &&
                        k_testdata[i][1] == 0)
                        continue;

                    acb_zero(z);
                    arf_set_d(arb_midref(acb_realref(z)), k_testdata[i][0]);
                    arf_set_d(arb_midref(acb_imagref(z)), cnj ? -k_testdata[i][1] : k_testdata[i][1]);

                    acb_modular_elliptic_k_cpx(w1, z, NUM_DERIVS, prec);

                    for (j = 0; j < NUM_DERIVS; j++)
                    {
                        arf_set_d(arb_midref(acb_realref(t)), k_testdata[i][2+2*j]);
                        mag_set_d(arb_radref(acb_realref(t)), fabs(k_testdata[i][2+2*j]) * EPS);
                        arf_set_d(arb_midref(acb_imagref(t)), cnj ? -k_testdata[i][2+2*j+1] : k_testdata[i][2+2*j+1]);
                        mag_set_d(arb_radref(acb_imagref(t)), fabs(k_testdata[i][2+2*j+1]) * EPS);

                        if (!acb_overlaps(w1 + j, t))
                        {
                            flint_printf("FAIL\n\n");
                            flint_printf("j = %wd\n\n", j);
                            flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
                            flint_printf("t = "); acb_printd(t, 15); flint_printf("\n\n");
                            flint_printf("w1 = "); acb_printd(w1 + j, 15); flint_printf("\n\n");
                            flint_abort();
                        }
                    }
                }
            }
        }

        _acb_vec_clear(w1, NUM_DERIVS);
        acb_clear(z);
        acb_clear(t);
    }

    /* self-consistency test */
    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr m1, m2;
        acb_t z1, z2, t;
        slong i, len1, len2, prec1, prec2;

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);

        m1 = _acb_vec_init(len1);
        m2 = _acb_vec_init(len2);

        acb_init(z1);
        acb_init(z2);
        acb_init(t);

        acb_randtest(z1, state, prec1, 1 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            acb_set(z2, z1);
        }
        else
        {
            acb_randtest(t, state, prec2, 1 + n_randint(state, 100));
            acb_add(z2, z1, t, prec2);
            acb_sub(z2, z2, t, prec2);
        }

        acb_modular_elliptic_k_cpx(m1, z1, len1, prec1);
        acb_modular_elliptic_k_cpx(m2, z2, len2, prec2);

        for (i = 0; i < FLINT_MIN(len1, len2); i++)
        {
            if (!acb_overlaps(m1 + i, m2 + i))
            {
                flint_printf("FAIL (overlap)\n\n");
                flint_printf("iter = %wd, i = %wd, len1 = %wd, len2 = %wd, prec1 = %wd, prec2 = %wd\n\n",
                    iter, i, len1, len2, prec1, prec2);

                flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n\n");
                flint_printf("m1 = "); acb_printd(m1, 30); flint_printf("\n\n");
                flint_printf("m2 = "); acb_printd(m2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        _acb_vec_clear(m1, len1);
        _acb_vec_clear(m2, len2);

        acb_clear(z1);
        acb_clear(z2);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
#undef NUM_TESTS
#undef EPS
