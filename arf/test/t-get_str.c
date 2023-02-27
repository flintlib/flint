/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "arf.h"
#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_str....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arf_t x;
        arb_t y;
        char *s1, * s2;
        slong d;

        arf_init(x);
        arb_init(y);

        arf_randtest_special(x, state, 2 + n_randint(state, 200), 2 + n_randint(state, 100));
        arb_set_arf(y, x);
        d = 1 + n_randint(state, 200);

        s1 = arf_get_str(x, d);
        s2 = arb_get_str(y, d, ARB_STR_NO_RADIUS);

        if ((arf_is_pos_inf(x) && strcmp(s1, "+inf")) ||
            (arf_is_neg_inf(x) && strcmp(s1, "-inf")) ||
            (arf_is_nan(x) && strcmp(s1, "nan")) ||
            (arf_is_finite(x) && strcmp(s1, s2)))
        {
            flint_printf("FAIL\n");
            arf_print(x); flint_printf("\n");
            arb_print(y); flint_printf("\n");
            flint_printf("%s\n", s1);
            flint_printf("%s\n", s2);
            flint_abort();
        }

        flint_free(s1);
        flint_free(s2);

        arf_clear(x);
        arb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
