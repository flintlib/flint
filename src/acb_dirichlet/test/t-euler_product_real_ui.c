/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

const signed char chi[8][6] = {
  {1, 1},
  {2, 0, 1},
  {3, 0, 1, 1},
  {3, 0, 1, -1},
  {4, 0, 1, 0, 1},
  {4, 0, 1, 0, -1},
  {5, 0, 1, 1, 1, 1},
  {5, 0, 1, -1, -1, 1},
};

const double L10[8] = {
  1.0009945751278180853371459589,
  1.0000170413630448254881839023,
  1.00097762319679253337038382667,
  0.999024291448866695201196896346,
  1.0000170413630448254881839023,
  0.999983164026196877405540729958,
  1.00099447262597359224857402038,
  0.999007468458940084215357132419
};

TEST_FUNCTION_START(acb_dirichlet_euler_product_real_ui, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t res1, res2;
        ulong s;
        slong prec1, prec2, accuracy;
        int choice, reciprocal1, reciprocal2;

        if (iter % 10 == 0)
        {
            s = n_randtest(state);
            prec1 = 2 + n_randint(state, 300);
            prec2 = 2 + n_randint(state, 300);
        }
        else
        {
            s = 6 + n_randint(state, 1 << n_randint(state, 12));
            prec1 = 2 + n_randint(state, 12 * s);
            prec2 = 2 + n_randint(state, 12 * s);
        }

        if (n_randint(state, 30)  == 0)
            prec1 = 2 + n_randint(state, 4000);

        choice = n_randint(state, 7);
        reciprocal1 = n_randint(state, 2);
        reciprocal2 = n_randint(state, 2);

        arb_init(res1);
        arb_init(res2);
        arb_randtest(res1, state, 200, 100);

        _acb_dirichlet_euler_product_real_ui(res1, s, chi[choice] + 1,
            chi[choice][0], reciprocal1, prec1);
        _acb_dirichlet_euler_product_real_ui(res2, s, chi[choice] + 1,
            chi[choice][0], reciprocal2, prec2);

        if (reciprocal1 != reciprocal2)
            arb_inv(res2, res2, prec2);

        if (!arb_overlaps(res1, res2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("s = %wu\n\n", s);
            flint_printf("chi: %d\n", choice);
            flint_printf("res1 = "); arb_printd(res1, prec1 / 3.33); flint_printf("\n\n");
            flint_printf("res2 = "); arb_printd(res2, prec2 / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        if (s >= 6 && prec1 < 2 * s * log(s))
        {
            accuracy = arb_rel_accuracy_bits(res1);

            if (accuracy < prec1 - 4)
            {
                flint_printf("FAIL: accuracy = %wd, prec = %wd\n\n", accuracy, prec1);
                flint_printf("res1 = "); arb_printd(res1, prec1 / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (s == 10)
        {
            arf_set_d(arb_midref(res2), L10[choice]);
            mag_set_d(arb_radref(res2), 1e-15);
            if (reciprocal1)
                arb_inv(res2, res2, 53);

            if (!arb_overlaps(res1, res2))
            {
                flint_printf("FAIL: overlap (2)\n\n");
                flint_printf("s = %wu\n\n", s);
                flint_printf("chi: %d\n", choice);
                flint_printf("res1 = "); arb_printd(res1, prec1 / 3.33); flint_printf("\n\n");
                flint_printf("res2 = "); arb_printd(res2, prec2 / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(res1);
        arb_clear(res2);
    }

    TEST_FUNCTION_END(state);
}
