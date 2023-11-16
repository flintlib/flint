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

static const double testdata_rj[16][10] = {
    {0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.77688623778582332014, 0.0},
    {2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0, 0.0, 0.14297579667156753833, 0.0},
    {2.0, 0.0, 3.0, 0.0, 4.0, 0.0, -1.0, 1.0, 0.13613945827770535204, -0.3820756162442716425},
    {0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 1.6490011662710884518, 0.0},
    {-1.0, 1.0, -1.0, -1.0, 1.0, 0.0, 2.0, 0.0, 0.94148358841220238083, 0.0},
    {0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.8260115229009316249, 1.22906619086434715},
    {-1.0, 1.0, -1.0, -1.0, 1.0, 0.0, -3.0, 1.0, -0.61127970812028172124, -1.068403839000680788},
    {-1.0, 1.0, -2.0, -1.0, 0.0, -1.0, -1.0, 1.0, 1.8249027393703805305, -1.2218475784827035855},
    /* additional tests where Carlson's algorithm may not be valid */
    {-1.0, -0.5, -10.0, -6.0, -10.0, -3.0, -5.0, 10.0, 0.128470516743927699, 0.102175950778504625},
    {1.987, 0.0, 4.463, -1.614, 0.0, 0.0, -3.965, 0.0, -0.341575118513811305, -0.394703757004268486},
    {0.3068, 0.0, -4.037, 0.632, 1.654, 0.0, -0.9609, 0.0, -1.14735199581485639, -0.134450158867472264},
    {0.3068, 0.0, -4.037, -0.632, 1.654, 0.0, -0.9609, 0.0, 1.758765901861727, -0.161002343366626892},
    {0.3068, 0.0, -4.037, 0.0632, 1.654, 0.0, -0.9609, 0.0, -1.17157627949475577, -0.069182614173988811},
    {0.3068, 0.0, -4.037, -0.0632, 1.654, 0.0, -0.9609, 0.0, 1.77940452391261626, 0.0388711305592447234},
    {0.3068, 0.0, -4.037, 0.00632, 1.654, 0.0, -0.9609, 0.0, -1.17337595670549633, -0.0623069224526925},
    {0.3068, 0.0, -4.037, -0.00632, 1.654, 0.0, -0.9609, 0.0, 1.77806722756403055, 0.0592749824572262329},
};

static const double testdata_rd[6][8] = {
    {0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 1.7972103521033883112, 0.0},
    {2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 0.16510527294261053349, 0.0},
    {0.0, 1.0, 0.0, -1.0, 2.0, 0.0, 0.65933854154219768919, 0.0},
    {0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.2708196271909686299, 2.7811120159520578776},
    {0.0, 0.0, -1.0, 1.0, 0.0, 1.0, -1.8577235439239060056, -0.96193450888838559989},
    {-2.0, -1.0, 0.0, -1.0, -1.0, 1.0, 1.8249027393703805305, -1.2218475784827035855},
};

TEST_FUNCTION_START(acb_elliptic_rj, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y, z, p, r1, r2;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        acb_init(x);
        acb_init(y);
        acb_init(z);
        acb_init(p);
        acb_init(r1);
        acb_init(r2);

        if (iter == 0)
        {
            slong k;

            for (k = 0; k < 16; k++)
            {
                acb_set_d_d(x, testdata_rj[k][0], testdata_rj[k][1]);
                acb_set_d_d(y, testdata_rj[k][2], testdata_rj[k][3]);
                acb_set_d_d(z, testdata_rj[k][4], testdata_rj[k][5]);
                acb_set_d_d(p, testdata_rj[k][6], testdata_rj[k][7]);
                acb_set_d_d(r2, testdata_rj[k][8], testdata_rj[k][9]);
                mag_set_d(arb_radref(acb_realref(r2)), 1e-14 * fabs(testdata_rj[k][8]));
                mag_set_d(arb_radref(acb_imagref(r2)), 1e-14 * fabs(testdata_rj[k][9]));

                for (prec1 = 16; prec1 <= 256; prec1 *= 2)
                {
                    acb_elliptic_rj(r1, x, y, z, p, 0, prec1);

                    if (!acb_overlaps(r1, r2) || (k < 8 && acb_rel_accuracy_bits(r1) < prec1 - 12))
                    {
                        flint_printf("FAIL: overlap (testdata rj)\n\n");
                        flint_printf("prec = %wd, accuracy = %wd\n\n", prec1, acb_rel_accuracy_bits(r1));
                        flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
                        flint_printf("y = "); acb_printd(y, 30); flint_printf("\n\n");
                        flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                        flint_printf("p = "); acb_printd(p, 30); flint_printf("\n\n");
                        flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                        flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                        flint_abort();
                    }
                }
            }

            for (k = 0; k < 6; k++)
            {
                acb_set_d_d(x, testdata_rd[k][0], testdata_rd[k][1]);
                acb_set_d_d(y, testdata_rd[k][2], testdata_rd[k][3]);
                acb_set_d_d(z, testdata_rd[k][4], testdata_rd[k][5]);
                acb_set_d_d(r2, testdata_rd[k][6], testdata_rd[k][7]);
                mag_set_d(arb_radref(acb_realref(r2)), 1e-14 * fabs(testdata_rd[k][6]));
                mag_set_d(arb_radref(acb_imagref(r2)), 1e-14 * fabs(testdata_rd[k][7]));

                for (prec1 = 16; prec1 <= 256; prec1 *= 2)
                {
                    acb_elliptic_rj(r1, x, y, z, z, 0, prec1);

                    if (!acb_overlaps(r1, r2) || acb_rel_accuracy_bits(r1) < prec1 - 12)
                    {
                        flint_printf("FAIL: overlap (testdata rd)\n\n");
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
        acb_set(p, z);

        /* Potentially slow general tests */
        if ((arb_is_nonnegative(acb_realref(x)) &&
              arb_is_nonnegative(acb_realref(y)) &&
              arb_is_nonnegative(acb_realref(z)) &&
              arb_is_positive(acb_realref(p))) || (n_randint(state, 10) == 0 && prec1 < 100 && prec2 < 100))
        {
            acb_elliptic_rj(r1, x, y, z, z, 0, prec1);
            acb_elliptic_rj(r2, x, y, z, p, 0, prec2);

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap R_D\n\n");
                flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
                flint_printf("y = "); acb_printd(y, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("p = "); acb_printd(p, 30); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }

            acb_randtest(p, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));

            /* Specialize */
            if (n_randint(state, 2))
            {
                if (arb_is_zero(acb_imagref(y)))
                    arb_one(acb_imagref(y));
                acb_conj(z, y);
                arb_zero(acb_imagref(x));
                arb_abs(acb_realref(x), acb_realref(x));
            }

            acb_elliptic_rj(r1, x, y, z, p, 0, prec1);

            switch (n_randint(state, 6))
            {
                case 0:
                    acb_elliptic_rj(r2, x, y, z, p, 0, prec2);
                    break;
                case 1:
                    acb_elliptic_rj(r2, x, z, y, p, 0, prec2);
                    break;
                case 2:
                    acb_elliptic_rj(r2, y, x, z, p, 0, prec2);
                    break;
                case 3:
                    acb_elliptic_rj(r2, y, z, x, p, 0, prec2);
                    break;
                case 4:
                    acb_elliptic_rj(r2, z, x, y, p, 0, prec2);
                    break;
                default:
                    acb_elliptic_rj(r2, z, y, x, p, 0, prec2);
                    break;
            }

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
                flint_printf("y = "); acb_printd(y, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("p = "); acb_printd(p, 30); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
        acb_clear(p);
        acb_clear(r1);
        acb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
