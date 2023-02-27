/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("asec_pi....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        slong p, p2;
        ulong q, q2;
        int res;

        qqbar_init(x);
        qqbar_init(y);

        do
        {
            q = 1 + n_randint(state, 15);
            p = n_randint(state, 1000);
            p -= 500;

            /* Don't generate poles */
        } while ((q / n_gcd(FLINT_ABS(p), q) == 2));

        qqbar_sec_pi(x, p, q);
        res = qqbar_asec_pi(&p2, &q2, x);
        if (res)
            qqbar_sec_pi(y, p2, q2);

        if (!res || !qqbar_equal(x, y) || n_gcd(FLINT_ABS(p2), q2) != 1 || !(0 <= p2 && p2 <= (slong) q2))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("res = %d\n\n", res);
            flint_printf("p, p2 = %wd %wd\n\n", p, p2);
            flint_printf("q, q2 = %wu %wu\n\n", q, q2);
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

