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

TEST_FUNCTION_START(acb_theta_char_is_syzygous, state)
{
    slong iter;

    /* Test: there are 60 syzygous triples for g = 2 */
    for (iter = 0; iter < 1; iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        ulong ch1, ch2, ch3;
        slong cnt = 0;

        for (ch1 = 0; ch1 < n; ch1++)
        {
            for (ch2 = ch1; ch2 < n; ch2++)
            {
                for (ch3 = ch2; ch3 < n; ch3++)
                {
                    if (acb_theta_char_is_syzygous(ch1, ch2, ch3, g))
                    {
                        cnt++;
                    }
                }
            }
        }

        if (cnt != 60)
        {
            flint_printf("FAIL (cnt = %wd)\n", cnt);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
