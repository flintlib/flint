/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly_factor.h"

/* manual slow implementation of zassenhaus_subset_next_disjoint */
static int my_next_disjoint(slong * s, slong r)
{
    slong i, j;
    slong * olds;

    olds = flint_malloc(r*sizeof(slong));

    for (i = 0; i < r; i++)
        olds[i] = s[i];

next:

    if (!zassenhaus_subset_next(s, r))
    {
        flint_free(olds);
        return 0;
    }

    for (i = 0; i < r; i++)
    {
        if (olds[i] >= 0 && s[i] >= 0)
            goto next;
    }

    j = 0;
    for (i = 0; i < r; i++)
    {
        if (olds[i] < 0)
            s[j++] = s[i];
    }

    flint_free(olds);

    return 1;
}

/* length of s is r */
static void my_subset_print(slong * s, slong r)
{
    slong i;

    FLINT_ASSERT(r >= 0);

    for (i = 0; i < r; i++)
    {
        if (s[i] >= 0)
            flint_printf(" %wd", s[i]);
        else
            flint_printf(" X", s[i]);
    }
}

TEST_FUNCTION_START(fmpz_poly_factor_zassenhaus_subset, state)
{
    fmpz_t f;
    slong i, j, k, r;
    slong * s, * s1, * s2;
    int res1, res2;

    fmpz_init(f);

    for (r = 0; r <= 13; r++)
    {
        s = (slong *) flint_malloc(r*sizeof(slong));
        s1 = (slong *) flint_malloc(r*sizeof(slong));
        s2 = (slong *) flint_malloc(r*sizeof(slong));

        for (j = 0; j < r; j++)
            s[j] = j;

        for (i = 0; i <= r; i++)
        {
            ulong cnt = 0;
            zassenhaus_subset_first(s, r, i);
            do {
                cnt++;
                k = 0;
                for (j = 0; j < r; j++)
                {
                    if (s[j] != j && s[j] != -j - 1)
                    {
                        flint_printf("FAIL\ncheck subset element\n");
                        fflush(stdout);
                        flint_abort();
                    }
                    k += (s[j] >= 0);
                }
                if (k != i)
                {
                    flint_printf("FAIL\ncheck subset size\n");
                    fflush(stdout);
                    flint_abort();
                }

                for (j = 0; j < r; j++)
                {
                    s1[j] = s[j];
                    s2[j] = s[j];
                }

                res1 = zassenhaus_subset_next_disjoint(s1, r);
                res2 = my_next_disjoint(s2, r);
                if (res1 != res2)
                {
                    flint_printf("FAIL\ncheck next disjoint success\n");
                    flint_printf("s:");
                    my_subset_print(s, r);
                    flint_printf("\n");
                    flint_printf("res1: %d  res2: %d\n", res1, res2);
                    fflush(stdout);
                    flint_abort();
                }

                k = 0;
                for (j = 0; j < r; j++)
                {
                    if (s[j] < 0)
                    {
                        if (s1[k] != s[j] && s1[k] != -s[j] - 1)
                        {
                            flint_printf("FAIL\ncheck next disjoint is disjoint\n");
                            fflush(stdout);
                            flint_abort();
                        }
                        k++;
                    }
                }

                if (res1)
                {
                    for (j = 0; j < r - i; j++)
                    {
                        if (s1[j] != s2[j])
                        {
                            flint_printf("FAIL\ncheck next disjoint subset\n");
                            flint_printf("s:");
                            my_subset_print(s, r);
                            flint_printf("\n");
                            flint_printf("s1: ");
                            my_subset_print(s1, r - i);
                            flint_printf("\n");
                            flint_printf("s2: ");
                            my_subset_print(s2, r - i);
                            flint_printf("\n");
                            fflush(stdout);
                            flint_abort();
                        }
                    }
                }
            } while (zassenhaus_subset_next(s, r));

            fmpz_bin_uiui(f, r, i);
            if (!fmpz_equal_ui(f, cnt))
            {
                flint_printf("FAIL\ncheck subset size\n");
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(s2);
        flint_free(s1);
        flint_free(s);
    }

    fmpz_clear(f);

    TEST_FUNCTION_END(state);
}
