/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    slong k, i, j, N, length, nfields, bits;
    ulong * a, * b, *max, *max2;
    ulong max_length, max_fields;
    FLINT_TEST_INIT(state);

    flint_printf("max_degrees....");
    fflush(stdout);

    max_length = 100;
    max_fields = 20;

    a    = (ulong *) flint_malloc(max_length*max_fields*sizeof(ulong));
    b    = (ulong *) flint_malloc(max_length*max_fields*sizeof(ulong));
    max  = (ulong *) flint_malloc(max_fields*sizeof(ulong));
    max2 = (ulong *) flint_malloc(max_fields*sizeof(ulong));

    for (k = 0; k < 100 * flint_test_multiplier(); k++)
    {
        /*
            calculate the maximum by hand using FLINT_BITS
            then pack FLINT_BITS => bits and compare the output
        */
        for (bits = 8; bits <= FLINT_BITS; bits += 1)
        {
            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, max_fields) + 1;
            N = words_per_exp(nfields, bits);

            for (j = 0; j < nfields; j++)
                max[j] = 0;

            for (i = 0; i < nfields*length; i += nfields)
                for (j = 0; j < nfields; j++)
                {
                    a[i + j] = n_randint(state, 0) & (l_shift(UWORD(1), bits)
                                                                          - 1);
                    max[nfields - j - 1] = FLINT_MAX(max[nfields - j - 1],
                                                                     a[i + j]);
                }

            /* FLINT_BITS => bits */
            for (i = 0; i < length; i++)
                mpoly_set_monomial(b + i*N, a + i*nfields, bits, nfields, 0, 0);

            mpoly_max_degrees(max2, b, length, bits, nfields);

            for (i = 0; i < nfields; i++)
                if (max[i] != max2[i])
                {
                    printf("FAIL\n");
                    flint_printf("bits = %wd, nfields = %wd\n", bits, nfields);
                    flint_abort();
                }
        }
    }

    flint_free(max2);
    flint_free(max);
    flint_free(b);
    flint_free(a);

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

