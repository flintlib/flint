/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("char_is_even....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: values for g = 2 */
    for (iter = 0; iter < 1; iter++)
    {
        slong g = 2;
        slong even[10] = {0, 1, 2, 3, 4, 6, 8, 9, 12, 15};
        slong odd[6] = {5, 7, 10, 11, 13, 14};
        slong i;

        for (i = 0; i < 10; i++)
        {
            if (!acb_theta_char_is_even(even[i], g))
            {
                flint_printf("FAIL (even)\n");
                flint_printf("i = %wd, ab = %wd\n", i, even[i]);
                flint_abort();
            }
        }

        for (i = 0; i < 6; i++)
        {
            if (acb_theta_char_is_even(odd[i], g))
            {
                flint_printf("FAIL (odd)\n");
                flint_printf("i = %wd, ab = %wd\n", i, odd[i]);
                flint_abort();
            }
        }
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
