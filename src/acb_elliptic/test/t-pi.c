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
#include "acb_modular.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

static const double testdata_pi[17][6] = {
  {-2.0, 0.0, 0.0, 0.0, 0.9068996821171089253, 0.0},
  {-2.0, 0.0, 2.0, 0.0, 0.90023699434287315449, -0.5446596026314193538},
  {-2.0, 0.0, 0.0, 1.0, 0.84791174714919877248, 0.12945491185218067805},
  {0.0, 0.0, -2.0, 0.0, 1.1714200841467698589, 0.0},
  {0.0, 0.0, 2.0, 0.0, 1.3110287771460599052, -1.3110287771460599052},
  {0.0, 0.0, 0.0, 1.0, 1.4212722810450360172, 0.29538028421477684284},
  {2.0, 0.0, -2.0, 0.0, 0.30328372333566606144, -1.1107207345395915618},
  {2.0, 0.0, 0.0, 0.0, 0.0, -1.5707963267948966192},
  {2.0, 0.0, 0.0, 1.0, 0.53486323549280133642, -1.7002121907485779385},
  {0.0, 1.0, -2.0, 0.0, 0.95136015191114896677, 0.33513918759114703301},
  {0.0, 1.0, 0.0, 0.0, 1.2203312255379458475, 0.50547774420519745381},
  {0.0, 1.0, 2.0, 0.0, 1.7987619218751253398, -0.55552744414075242238},
  {2.0, -1.0, 2.0, 1.0, 1.8578723151271115, -1.18642180609983531},
  {2.0, -0.5, 0.5, 1.0, 0.936761970766645807, -1.61876787838890786},
  {2.0, 0.0, 1.0, 1.0, 0.999881420735506708, -2.4139272867045391},
  {2.0, 1.0, 2.0, -1.0, 1.8578723151271115, 1.18642180609983531},
  {2.0, 1.0, 2.0, 0.0, 2.78474654927885845, 2.02204728966993314},
};

TEST_FUNCTION_START(acb_elliptic_pi, state)
{
    slong iter;

    /* Test self-consistency, and Pi(n,n) = E(n) / (1-n) */
    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t n, m, r1, r2, t;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        acb_init(n);
        acb_init(m);
        acb_init(r1);
        acb_init(r2);
        acb_init(t);

        if (iter == 0)
        {
            slong k;

            for (k = 0; k < 17; k++)
            {
                acb_set_d_d(n, testdata_pi[k][0], testdata_pi[k][1]);
                acb_set_d_d(m, testdata_pi[k][2], testdata_pi[k][3]);
                acb_set_d_d(r2, testdata_pi[k][4], testdata_pi[k][5]);
                mag_set_d(arb_radref(acb_realref(r2)), 1e-14 * fabs(testdata_pi[k][4]));
                mag_set_d(arb_radref(acb_imagref(r2)), 1e-14 * fabs(testdata_pi[k][5]));

                for (prec1 = 32; prec1 <= (0.1 * flint_test_multiplier() < 1.0 ? 64 : 256); prec1 *= 2)
                {
                    acb_elliptic_pi(r1, n, m, prec1 + 30);

                    if (!acb_overlaps(r1, r2) || acb_rel_accuracy_bits(r1) < prec1)
                    {
                        flint_printf("FAIL: overlap (testdata)\n\n");
                        flint_printf("prec = %wd, accuracy = %wd\n\n", prec1, acb_rel_accuracy_bits(r1));
                        flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
                        flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
                        flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                        flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                        flint_abort();
                    }
                }
            }
        }

        acb_randtest(n, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));
        acb_randtest(m, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));

        acb_one(t);
        if ((acb_is_real(n) && acb_is_real(m) &&
            arb_le(acb_realref(n), acb_realref(t)) &&
            arb_le(acb_realref(m), acb_realref(t))) || n_randint(state, 10) == 0)
        {
            acb_elliptic_pi(r1, n, n, prec1);

            acb_modular_elliptic_e(r2, n, prec1);
            acb_sub_ui(t, n, 1, prec1);
            acb_neg(t, t);
            acb_div(r2, r2, t, prec1);

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap (m=n)\n\n");
                flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }

            acb_elliptic_pi(r1, n, m, prec1);

            acb_randtest(t, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));
            acb_add(n, n, t, prec2);
            acb_sub(n, n, t, prec2);
            acb_randtest(t, state, 1 + n_randint(state, 300), 1 + n_randint(state, 30));
            acb_add(m, m, t, prec2);
            acb_sub(m, m, t, prec2);

            acb_elliptic_pi(r2, n, m, prec2);

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap (consistency)\n\n");
                flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
                flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(n);
        acb_clear(m);
        acb_clear(r1);
        acb_clear(r2);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
