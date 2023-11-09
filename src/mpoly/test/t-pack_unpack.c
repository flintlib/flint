/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "mpoly.h"

TEST_FUNCTION_START(mpoly_pack_unpack, state)
{
    slong k, i, length, nfields, bits1, bits2;
    ulong * a, * b, * c, * d;
    ulong max_length, max_fields;

    max_length = 50;
    max_fields = 20;

    a = flint_malloc(max_length*max_fields*sizeof(ulong));
    b = flint_malloc(max_length*max_fields*sizeof(ulong));
    c = flint_malloc(max_length*max_fields*sizeof(ulong));
    d = flint_malloc(max_length*max_fields*sizeof(ulong));

    for (k = 0; k < 10 * flint_test_multiplier(); k++)
    {
        /* do FLINT_BITS => bits1 => FLINT_BITS and compare */
        for (bits1 = 8; bits1 <= FLINT_BITS; bits1 += 1)
        {
            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, max_fields) + 1;

            for (i = 0; i < length*nfields; i++)
                a[i] = n_randint(state, 0) & (l_shift(UWORD(1), bits1) - 1);

            mpoly_pack_vec_ui(b, a, bits1, nfields, length);
            mpoly_unpack_vec_ui(c, b, bits1, nfields, length);

            for (i = 0; i < length*nfields; i++)
                if (a[i] != c[i])
                {
                    printf("FAIL\n");
                    flint_printf("bits1 = %wd  fields = %wd  len = %wd\n", bits1, nfields, length);
                    fflush(stdout);
                    flint_abort();
                }
        }
    }

    for (k = 0; k < 1 * flint_test_multiplier(); k++)
    {
        /* do FLINT_BITS => bits1 => bits2 => FLINT_BITS and compare */
        for (bits1 = MPOLY_MIN_BITS; bits1 <= FLINT_BITS; bits1 += 1)
        {
        for (bits2 = MPOLY_MIN_BITS; bits2 <= FLINT_BITS; bits2 += 1)
        {
            ulong mask = (l_shift(UWORD(1), bits1 - 1) - 1)
                       & (l_shift(UWORD(1), bits2 - 1) - 1);
            mpoly_ctx_t mctx;

            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, max_fields) + 1;

            mpoly_ctx_init(mctx, nfields, ORD_LEX);

            for (i = 0; i < length*nfields; i++)
            {
                /* leave room for sign bit in the repacking */
                a[i] = n_randlimb(state) & mask;
            }

            mpoly_repack_monomials(b, bits1,      a, FLINT_BITS, length, mctx);
            mpoly_repack_monomials(c, bits2,      b, bits1,      length, mctx);
            mpoly_repack_monomials(d, FLINT_BITS, c, bits2,      length, mctx);

            for (i = 0; i < length*nfields; i++)
            {
                if (a[i] != d[i])
                {
                    printf("FAIL\n");
                    flint_printf("bits1 = %wd, bits2 = %wd\n", bits1, bits2);
                    fflush(stdout);
                    flint_abort();
                }
            }

            mpoly_ctx_clear(mctx);
        }
        }
    }

    flint_free(d);
    flint_free(c);
    flint_free(b);
    flint_free(a);

    TEST_FUNCTION_END(state);
}
