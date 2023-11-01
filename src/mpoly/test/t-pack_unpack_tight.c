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

TEST_FUNCTION_START(mpoly_pack_unpack_tight, state)
{
    slong k, i, j, length, nfields, bits1, bits2;
    slong * bases;
    ulong * a, * b, * c, * d, * t;
    ulong max_length, max_fields;

    max_length = 50;
    max_fields = FLINT_BITS/8;  /* exponents should fit in one word */

    a = flint_malloc(max_length*max_fields*sizeof(ulong));
    b = flint_malloc(max_length*max_fields*sizeof(ulong));
    c = flint_malloc(max_length*max_fields*sizeof(ulong));
    d = flint_malloc(max_length*max_fields*sizeof(ulong));
    t = flint_malloc(max_length*sizeof(ulong));
    bases = flint_malloc(max_fields*sizeof(slong));

    for (k = 0; k < 20 * flint_test_multiplier(); k++)
    {
        /* do FLINT_BITS => bits1
                         => tight packing
                         => bits2 => FLINT_BITS and compare */

        for (bits1 = 8; bits1 <= FLINT_BITS; bits1 += 1)
        for (bits2 = 8; bits2 <= FLINT_BITS; bits2 += 1)
        {
            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, FLINT_BITS/FLINT_MAX(bits1, bits2)) + 1;

            for (j = 0; j < nfields; j++)
                bases[j] =  n_randint(state, 200) + 1;

            for (i = 0; i < length; i += 1)
                for (j = 0; j < nfields; j++)
                    a[nfields*i + j] = n_randint(state, bases[j]);

            mpoly_pack_vec_ui(b, a, bits1, nfields, length);
            mpoly_pack_monomials_tight(t, b, length, bases, nfields, bits1);
            mpoly_unpack_monomials_tight(c, t, length, bases, nfields, bits2);
            mpoly_unpack_vec_ui(d, c, bits2, nfields, length);

            for (i = 0; i < length*nfields; i++)
                if (a[i] != d[i])
                {
                    printf("FAIL\nunpack_monomials_tight\n");
                    flint_printf("bits1 = %wd, bits2 = %wd\n", bits1, bits2);
                    fflush(stdout);
                    flint_abort();
                }
        }
    }

    flint_free(bases);
    flint_free(t);
    flint_free(d);
    flint_free(c);
    flint_free(b);
    flint_free(a);

    TEST_FUNCTION_END(state);
}
