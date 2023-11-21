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

TEST_FUNCTION_START(acb_theta_char_is_goepel, state)
{
    slong iter;

    /* Test: there are 15 Goepel quadruples for g = 2 */
    for (iter = 0; iter < 1; iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        ulong ch1, ch2, ch3, ch4;
        slong cnt = 0;

        for (ch1 = 0; ch1 < n; ch1++)
        {
            for (ch2 = ch1; ch2 < n; ch2++)
            {
                for (ch3 = ch2; ch3 < n; ch3++)
                {
                    for (ch4 = ch3; ch4 < n; ch4++)
                    {
                        if (acb_theta_char_is_goepel(ch1, ch2, ch3, ch4, g))
                        {
                            cnt++;
                        }
                    }
                }
            }
        }

        if (cnt != 15)
        {
            flint_printf("FAIL (cnt = %wd)\n", cnt);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
