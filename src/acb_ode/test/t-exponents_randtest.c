/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Generated using Claude Opus 4.8 */

#include "test_helpers.h"
#include "acb.h"
#include "acb_ode.h"

/* tests exponents_randtest + exponents accessors and group_set */

static void
check_group(const acb_ode_group_t grp, slong disp)
{
    slong glen = acb_ode_group_length(grp);

    FLINT_TEST(acb_ode_group_nlogs(grp, WORD_MAX) == glen);
    FLINT_TEST(acb_ode_group_nlogs(grp, disp) == glen);
    FLINT_TEST(acb_ode_group_nlogs(grp, -1) == 0);
    FLINT_TEST(acb_ode_group_multiplicity(grp, -1) == 0);

    slong prev = 0, cumul = 0;
    for (slong n = 0; n <= disp; n++)
    {
        slong cur = acb_ode_group_nlogs(grp, n);
        slong mult = acb_ode_group_multiplicity(grp, n);

        FLINT_TEST(cur >= prev);
        FLINT_TEST(cur - prev == mult);

        cumul += mult;
        prev = cur;
    }
    FLINT_TEST(cumul == glen);
}

TEST_FUNCTION_START(acb_ode_exponents_randtest, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        acb_ode_exponents_t expos;
        slong len, disp, prec, bits;

        acb_ode_exponents_init(expos);

        len = n_randint(state, 10);
        disp = n_randint(state, 30);
        prec = 2 + n_randint(state, 80);
        bits = 1 + n_randint(state, 10);

        if (n_randint(state, 2))
            /* test overwriting an existing object */
            acb_ode_exponents_randtest(expos, state, 1 + n_randint(state, 6),
                                       10, prec, 3);
        acb_ode_exponents_randtest(expos, state, len, disp, prec, bits);

        FLINT_TEST(acb_ode_exponents_length(expos) == len);

        slong total = 0;
        for (slong g = 0; g < expos->ngroups; g++)
        {
            check_group(expos->groups + g, disp);
            total += acb_ode_group_length(expos->groups + g);
        }
        FLINT_TEST(total == len);

        if (expos->ngroups > 0)
        {
            slong g = n_randint(state, expos->ngroups);
            acb_ode_group_struct * src = expos->groups + g;
            acb_ode_group_t dst;

            acb_ode_group_init(dst, src->nshifts);
            acb_ode_group_set(dst, src);

            FLINT_TEST(acb_equal(dst->leader, src->leader));
            FLINT_TEST(acb_ode_group_length(dst) == acb_ode_group_length(src));
            for (slong n = -1; n <= disp; n++)
            {
                FLINT_TEST(acb_ode_group_multiplicity(dst, n)
                           == acb_ode_group_multiplicity(src, n));
                FLINT_TEST(acb_ode_group_nlogs(dst, n)
                           == acb_ode_group_nlogs(src, n));
            }

            acb_ode_group_clear(dst);
        }

        acb_ode_exponents_clear(expos);
    }

    TEST_FUNCTION_END(state);
}
