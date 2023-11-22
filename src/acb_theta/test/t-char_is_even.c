/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_char_is_even, state)
{
    slong iter;

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

    TEST_FUNCTION_END(state);
}
