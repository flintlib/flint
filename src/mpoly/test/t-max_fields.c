/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpoly.h"

TEST_FUNCTION_START(mpoly_max_fields, state)
{
    slong k, i, j, length, nfields, bits;
    ulong * a, * b, * max, * max2;
    ulong max_length, max_fields;

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
            mpoly_ctx_t mctx;
            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, max_fields) + 1;

            mpoly_ctx_init(mctx, nfields, ORD_LEX);

            for (j = 0; j < nfields; j++)
                max[j] = 0;

            for (i = 0; i < nfields*length; i += nfields)
            {
                for (j = 0; j < nfields; j++)
                {
                    a[i + j] = n_randint(state, 0);
                    a[i + j] &= (UWORD(1) << (bits - 1)) - 1;
                    max[j] = FLINT_MAX(max[j], a[i + j]);
                }
            }

            /* FLINT_BITS => bits */
            mpoly_pack_vec_ui(b, a, bits, nfields, length);

            mpoly_max_fields_ui_sp(max2, b, length, bits, mctx);

            for (i = 0; i < nfields; i++)
                if (max[i] != max2[i])
                {
                    printf("FAIL\n");
                    flint_printf("bits = %wd, nfields = %wd\n", bits, nfields);
                    fflush(stdout);
                    flint_abort();
                }
        }
    }

    flint_free(max2);
    flint_free(max);
    flint_free(b);
    flint_free(a);

    TEST_FUNCTION_END(state);
}
