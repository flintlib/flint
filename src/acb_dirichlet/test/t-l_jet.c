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
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

/* Laurent expansions at s = 1 of first 10 principal L-functions */

/* with mpmath:
chis = [[1],[0,1],[0,1,1],[0,1,0,1],[0,1,1,1,1],[0,1,0,0,0,1],[0,1,1,1,1,1,1],
 [0,1,0,1,0,1,0,1],[0,1,1,0,1,1,0,1,1],[0,1,0,1,0,0,0,1,0,1]]

mp.dps = 40
for chi in chis:
    phi = chi.count(1); q = len(chi)
    L = lambda s: dirichlet(s, chi) - phi/((s-1)*q)
    c0 = taylor(L, 1, 0, method="quad")
    c1 = taylor(L, 1, 5, singular=True)[1:]
    for c in c0 + c1:
        print nstr(c, 20) + ",",
    print
*/

#define TESTQ 10
#define TESTLEN 6

static const double laurent_data[TESTQ][TESTLEN] = {
    {0.57721566490153286061, 0.072815845483676724861, -0.0048451815964361592423,
        -0.00034230573671722431103, 0.000096890419394470835728, -6.6110318108421891813e-6},
    {0.63518142273073908501, 0.11634237461305384831, -0.018765738937942729408,
        0.00061334298434914532242, 0.00042338142025747308027, -0.00010545096447379519004},
    {0.7510145394903918042, 0.058764477744540050414, -0.019011359100973296683,
        0.0056382252365739175151, -0.0009550480622176659462, 0.000021808301216554848718},
    {0.63518142273073908501, 0.11634237461305384831, -0.018765738937942729408,
        0.00061334298434914532242, 0.00042338142025747308027, -0.00010545096447379519004},
    {0.78366011440804636341, -0.014977808062405260803, 0.0090104707969118845102,
        0.003603799084856807634, -0.0029351216034181476022, 0.00093077685173004747355},
    {0.60655632993184433857, 0.2095885418562151802, -0.060844893711330538429,
        0.0068080382961291386117, 0.0022236616427578346453, -0.0013581825996235430782},
    {0.77274344835207292411, -0.047596894381510269689, 0.035406039531261788462,
        -0.0054159870134630085898, -0.0019749752308692423114, 0.0014492998471928196325},
    {0.63518142273073908501, 0.11634237461305384831, -0.018765738937942729408,
        0.00061334298434914532242, 0.00042338142025747308027, -0.00010545096447379519004},
    {0.7510145394903918042, 0.058764477744540050414, -0.019011359100973296683,
        0.0056382252365739175151, -0.0009550480622176659462, 0.000021808301216554848718},
    {0.66908892942800130547, 0.16801639259476784034, -0.072611999814034642781,
        0.024624650443138705595, -0.004951850872731033514, -0.00020178815459414925709}
};

TEST_FUNCTION_START(acb_dirichlet_l_jet, state)
{
    slong iter;

    /* test Laurent series at s = 1 */
    {
        acb_t s, t;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        acb_ptr vec;
        ulong q;
        slong i;

        acb_init(s);
        acb_init(t);
        vec = _acb_vec_init(TESTLEN);

        acb_one(s);
        for (q = 1; q <= TESTQ; q++)
        {
            dirichlet_group_init(G, q);
            dirichlet_char_init(chi, G);

            acb_dirichlet_l_jet(vec, s, G, chi, 1, TESTLEN, 100);

            for (i = 0; i < TESTLEN; i++)
            {
                acb_set_d(t, laurent_data[q - 1][i]);
                mag_set_d(arb_radref(acb_realref(t)),
                    fabs(laurent_data[q - 1][i]) * 1e-14);

                if (!acb_overlaps(vec + i, t))
                {
                    flint_printf("FAIL: Laurent series\n\n");
                    flint_printf("q = %wu  i = %wd\n\n", q, i);
                    flint_printf("r1 = "); acb_printn(vec + i, 50, 0); flint_printf("\n\n");
                    flint_printf("r2 = "); acb_printn(t, 50, 0); flint_printf("\n\n");
                    flint_abort();
                }
            }

            dirichlet_char_clear(chi);
            dirichlet_group_clear(G);
        }

        acb_clear(s);
        acb_clear(t);
        _acb_vec_clear(vec, TESTLEN);
    }

    /* test self-consistency */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        acb_ptr vec1, vec2;
        slong len1, len2;
        slong prec1, prec2;
        int deflate1, deflate2;
        ulong q, k;
        slong i;

        len1 = n_randint(state, 5);
        len2 = n_randint(state, 5);
        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);
        deflate1 = n_randint(state, 2);
        deflate2 = n_randint(state, 2);
        q = 1 + n_randint(state, 20);
        k = n_randint(state, n_euler_phi(q));

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, k);

        acb_init(s);
        vec1 = _acb_vec_init(len1);
        vec2 = _acb_vec_init(len2);

        if (n_randint(state, 4) == 0)
            acb_one(s);
        else
            acb_randtest(s, state, 2 + n_randint(state, 200), 2);

        acb_dirichlet_l_jet(vec1, s, G, chi, deflate1, len1, prec1);
        acb_dirichlet_l_jet(vec2, s, G, chi, deflate2, len2, prec2);

        if (deflate1 != deflate2 && dirichlet_char_is_principal(G, chi))
        {
            /* add or subtract phi(q)/((s+x-1)q) */
            acb_t t, u;
            acb_init(t);
            acb_init(u);

            acb_set_ui(t, n_euler_phi(q));
            acb_div_ui(t, t, q, prec1);
            acb_sub_ui(u, s, 1, prec1);

            for (i = 0; i < len1; i++)
            {
                acb_div(t, t, u, prec1);
                if (deflate1)
                    acb_add(vec1 + i, vec1 + i, t, prec1);
                else
                    acb_sub(vec1 + i, vec1 + i, t, prec1);
                acb_neg(t, t);
            }

            acb_clear(t);
            acb_clear(u);
        }

        for (i = 0; i < FLINT_MIN(len1, len2); i++)
        {
            if (!acb_overlaps(vec1 + i, vec2 + i))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd  q = %wu  k = %wu  i = %wd\n\n", iter, q, k, i);
                flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printn(vec1 + i, 50, 0); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printn(vec2 + i, 50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        acb_clear(s);
        _acb_vec_clear(vec1, len1);
        _acb_vec_clear(vec2, len2);
    }

    TEST_FUNCTION_END(state);
}
