/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_elliptic.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

/* Test input from Carlson's paper and checked with mpmath. */

static const double testdata_rg[7][8] = {
    {0.0, 0.0, 16.0, 0.0, 16.0, 0.0, 3.1415926535897932385, 0.0},
    {2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 1.7255030280692277601, 0.0},
    {0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.4236065423969895433, 0.0},
    {-1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.44660591677018372657, 0.70768352357515390073},
    {0.0, -1.0, -1.0, 1.0, 0.0, 1.0, 0.36023392184473309034, 0.40348623401722113741},
    {0.0, 0.0, 0.0796, 0.0, 4.0, 0.0, 1.0284758090288040022, 0.0},
    /* more tests */
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
};

TEST_FUNCTION_START(acb_elliptic_rg, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y, z, r1, r2;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        acb_init(x);
        acb_init(y);
        acb_init(z);
        acb_init(r1);
        acb_init(r2);

        if (iter == 0)
        {
            slong k;

            for (k = 0; k < 7; k++)
            {
                acb_set_d_d(x, testdata_rg[k][0], testdata_rg[k][1]);
                acb_set_d_d(y, testdata_rg[k][2], testdata_rg[k][3]);
                acb_set_d_d(z, testdata_rg[k][4], testdata_rg[k][5]);
                acb_set_d_d(r2, testdata_rg[k][6], testdata_rg[k][7]);
                mag_set_d(arb_radref(acb_realref(r2)), 1e-14 * fabs(testdata_rg[k][6]));
                mag_set_d(arb_radref(acb_imagref(r2)), 1e-14 * fabs(testdata_rg[k][7]));

                for (prec1 = 16; prec1 <= 256; prec1 *= 2)
                {
                    acb_elliptic_rg(r1, x, y, z, 0, prec1);

                    if (!acb_overlaps(r1, r2) || acb_rel_accuracy_bits(r1) < prec1 * 0.9 - 10)
                    {
                        flint_printf("FAIL: overlap (testdata rg)\n\n");
                        flint_printf("prec = %wd, accuracy = %wd\n\n", prec1, acb_rel_accuracy_bits(r1));
                        flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
                        flint_printf("y = "); acb_printd(y, 30); flint_printf("\n\n");
                        flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                        flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                        flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                        flint_abort();
                    }
                }
            }
        }

        acb_randtest(x, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));
        acb_randtest(y, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));
        acb_randtest(z, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));

        acb_elliptic_rg(r1, x, y, z, 0, prec1);

        switch (n_randint(state, 6))
        {
            case 0:
                acb_elliptic_rg(r2, x, y, z, 0, prec2);
                break;
            case 1:
                acb_elliptic_rg(r2, x, z, y, 0, prec2);
                break;
            case 2:
                acb_elliptic_rg(r2, y, x, z, 0, prec2);
                break;
            case 3:
                acb_elliptic_rg(r2, y, z, x, 0, prec2);
                break;
            case 4:
                acb_elliptic_rg(r2, z, x, y, 0, prec2);
                break;
            default:
                acb_elliptic_rg(r2, z, y, x, 0, prec2);
                break;
        }

        if (!acb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); acb_printd(y, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
        acb_clear(r1);
        acb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
