/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"
#include "math.h"

int main(void)
{
    fmpz_t hh, kk;
    fmpq_t s1, s2;
    slong h, k;

    FLINT_TEST_INIT(state);

    flint_printf("dedekind_sum_coprime_large....");
    fflush(stdout);

    fmpz_init(hh);
    fmpz_init(kk);
    fmpq_init(s1);
    fmpq_init(s2);

    for (k = 0; k < 300; k++)
    {
        for (h = 0; h <= k; h++)
        {
            if (n_gcd(k, h) == 1)
            {
                fmpz_set_si(hh, h);
                fmpz_set_si(kk, k);

                fmpq_dedekind_sum_coprime_large(s1, hh, kk);
                fmpq_dedekind_sum_naive(s2, hh, kk);

                if (!fmpq_equal(s1, s2))
                {
                    flint_printf("FAIL:\n");
                    flint_printf("s(%wd,%wd)\n", h, k);
                    flint_printf("s1: "); fmpq_print(s1); flint_printf("\n");
                    flint_printf("s2: "); fmpq_print(s2); flint_printf("\n");
                    abort();
                }
            }
        }
    }

    fmpz_clear(hh);
    fmpz_clear(kk);
    fmpq_clear(s1);
    fmpq_clear(s2);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
