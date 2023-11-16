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

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

#define EPS 1e-13
#define NUM_DERIVS 4
#define NUM_TESTS 11

const double agm1_testdata[NUM_TESTS][10] = {
    {1.0, 0.0, 1.0, 0.0, 0.5, 0.0, -0.0625, 0.0, 0.03125, 0.0},
    {0.25, 0.0, 0.56075714507190064253, 0.0, 0.76633907325304843764, 0.0,
        -0.58010113691169132987, 0.0, 1.2991960360521313649, 0.0},
    {1.0, 1.0, 1.0491605287327802205, 0.47815574608816122933,
        0.44643105633549979073, -0.08043578677710866283,
        -0.015455284495904882924, 0.031976374729700173479,
        -0.005073437378084728324, -0.010673958729796444985},
    {0.0, 1.0, 0.59907011736779610372, 0.59907011736779610372,
        0.43640657965245804105, -0.16266353771533806267,
        0.031271486774469549792, 0.031271486774469549792,
        -0.0084910439043492636266, 0.022780442870120286166},
    {-1.0, 1.0, 0.18841106798868002055, 0.77800407878895828015,
        0.39320630832295102335, -0.19287323123455182026,
        0.016808115488846724979, 0.020163502567351546742,
        -0.011189063384352573532, -0.0045629824424054816356},
    {-0.25, 0.0, 0.24392673474953340413, 0.27799893427725564501,
        0.079794586546549867566, -0.32574453868033984617,
        -0.27530931370867131683, -0.60287933825809104112,
        -0.54714813521966247214, -0.97291617788393560861},
    {-2.0, 0.0, -0.42296620840880168736, 0.66126618346180476447,
        0.29655367830470777795, -0.27314834694816402295,
        0.018331225229748304014, -0.043277191398267359935,
        0.0094104572902577354423, -0.021071819981331905571},
    {-0.99999994039535522461, 0.0, 0.0069953999943208591971,
        0.08334545077423930704, 12454.444757282471906, 73670.447089584396744,
        -87884330303.90316493, -553342797575.05204396, 909784741818095264.43,
        5883796037072491540.4},
    {-1.0000000596046447754, 0.0, -0.0069954003670321135601,
        0.083345455480285282001, 12454.451634795395449, -73670.529975671147878,
        87884324285.574775284, -553342761339.97964199, 909784713383817144.98,
        -5883795855289758433.3},
    {-1.0, 5.9604644775390625e-8, 2.3677293150757997928e-9,
        0.083932595219022127005, 75242.17179390509748, -0.041726661532440829443,
        -18488.306894081923308, 563725525035.34898339, -5988414396610684162.5,
        -92537341636.835656296},
    {-1.0, -5.9604644775390625e-8, 2.3677293150757997928e-9,
        -0.083932595219022127005, 75242.17179390509748, 0.041726661532440829443,
        -18488.306894081923308, -563725525035.34898339, -5988414396610684162.5,
        92537341636.835656296},
};

TEST_FUNCTION_START(acb_agm1, state)
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
                    if (cnj == 1 && agm1_testdata[i][0] < 0 &&
                        agm1_testdata[i][1] == 0)
                        continue;

                    acb_zero(z);
                    arf_set_d(arb_midref(acb_realref(z)), agm1_testdata[i][0]);
                    arf_set_d(arb_midref(acb_imagref(z)), cnj ? -agm1_testdata[i][1] : agm1_testdata[i][1]);

                    acb_agm1_cpx(w1, z, NUM_DERIVS, prec);

                    for (j = 0; j < NUM_DERIVS; j++)
                    {
                        arf_set_d(arb_midref(acb_realref(t)), agm1_testdata[i][2+2*j]);
                        mag_set_d(arb_radref(acb_realref(t)), fabs(agm1_testdata[i][2+2*j]) * EPS);
                        arf_set_d(arb_midref(acb_imagref(t)), cnj ? -agm1_testdata[i][2+2*j+1] : agm1_testdata[i][2+2*j+1]);
                        mag_set_d(arb_radref(acb_imagref(t)), fabs(agm1_testdata[i][2+2*j+1]) * EPS);

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
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
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

        acb_agm1_cpx(m1, z1, len1, prec1);
        acb_agm1_cpx(m2, z2, len2, prec2);

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
#undef EPS
#undef NUM_DERIVS
#undef NUM_TESTS
