/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qfb.h"

TEST_FUNCTION_START(qfb_reduced_forms, state)
{
    int result;
    qfb * forms;
    qfb * forms2;
    slong i, j, k, num, num2;

    for (i = 1; i < 30000; i++)
    {
        num = qfb_reduced_forms(&forms, -i);
        num2 = qfb_reduced_forms_large(&forms2, -i);

        result = (num == num2);
        for (j = 0; result == 1 && j < num; j++)
        {
            for (k = 0; k < num; k++)
               if (qfb_equal(forms + j, forms2 + k))
                   break;

            result &= (k != num);

            if (!result) break;
        }

        if (!result)
        {
            printf("FAIL:\n");

            if (num == num2)
            {
                qfb_print(forms + j);
                printf(" != ");
                qfb_print(forms2 + j);
                printf("\n");
             } else
                flint_printf("%wd != %wd\n", num, num2);

            flint_abort();
        }

        qfb_array_clear(&forms, num);
        qfb_array_clear(&forms2, num2);
    }

    TEST_FUNCTION_END(state);
}
