/*
    Copyright (C) 2017 Apoorv Mishra

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
main(void)
{
    ulong limit, rand_num;
    slong deviation;
    int i, j, result;

    FLINT_TEST_INIT(state);

    flint_printf("n_urandint....");
    fflush(stdout);

    /* Test for limit <= 1000 */
    for (limit = 1; limit <= 100 * flint_test_multiplier(); limit++)
    {
        int count[1000] = {0};

        for (i = 0; i < 1000 * flint_test_multiplier(); i++)
        {
            rand_num = n_urandint(state, limit);
            count[rand_num]++;
        }

        result = 1;
        for (i = 0; i < limit; i++)
        {
            deviation = count[i] - 10000/limit;
            if (deviation >= 315 || deviation <= -315)
            {
                result = 0;
                break;
            }
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("limit = %wu, deviation = %wd\n", limit, deviation);
            abort();
        }
    }

    /* Test for larger values of limit */
    for (i = 0; i < 1000; i++)
    {
        ulong count_in_subrange[4] = {0};

        limit = UWORD_MAX/(i + 2)*(i + 1);

        for (j = 0; j < 1000 * flint_test_multiplier(); j++)
        {
            rand_num = n_urandint(state, limit);

            if (rand_num >= 3*(limit >> 2))
            {
                count_in_subrange[3]++;
            }
            else if (rand_num >= 2*(limit >> 2))
            {
                count_in_subrange[2]++;
            }
            else if (rand_num >= (limit >> 2))
            {
                count_in_subrange[1]++;
            }
            else
            {
                count_in_subrange[0]++;
            }
        }

        result = 1;
        for (j = 0; j < 4; j++)
        {
            deviation = count_in_subrange[j] - (10000 >> 2);
            if (deviation >= WORD(575) || deviation <= WORD(-575))
            {
                result = 0;
                break;
            }
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("limit = %wu, deviation = %wd\n", limit, deviation);
            for (j = 0; j < 4; j++)
            {
                flint_printf("count_in_subrange[%d] = %wu\n", j, count_in_subrange[j]);
            }
            abort();
        }
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
