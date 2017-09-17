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
    slong k, i, N, length, nfields, bits1, bits2;
    ulong * a, * b, * c, * d;
    ulong max_length, max_fields;
    FLINT_TEST_INIT(state);

    flint_printf("pack_unpack....");
    fflush(stdout);

    max_length = 100;
    max_fields = 20;

    a = flint_malloc(max_length*max_fields*sizeof(ulong));
    b = flint_malloc(max_length*max_fields*sizeof(ulong));
    c = flint_malloc(max_length*max_fields*sizeof(ulong));
    d = flint_malloc(max_length*max_fields*sizeof(ulong));


    for (k = 0; k < 100 * flint_test_multiplier(); k++)
    {
        /* do FLINT_BITS => bits1 => FLINT_BITS and compare */
        for (bits1 = 8; bits1 <= FLINT_BITS; bits1 += 1)
        {
            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, max_fields) + 1;
            N = words_per_exp(nfields, bits1);
            for (i = 0; i < length*nfields; i++)
                a[i] = n_randint(state, 0) & (l_shift(UWORD(1), bits1) - 1);

            mpoly_pack_vec(b, a, bits1, nfields, length);
            mpoly_unpack_vec(c, b, bits1, nfields, length);

            for (i = 0; i < length*nfields; i++)
                if (a[i] != c[i])
                {
                    printf("FAIL\n");
                    flint_printf("bits1 = %wd  fields = %wd  len = %wd\n", bits1, nfields, length);
                    flint_abort();
                }
        }
    }


    for (k = 0; k < 10 * flint_test_multiplier(); k++)
    {
        /* do FLINT_BITS => bits1 => bits2 => FLINT_BITS and compare */
        for (bits1 = 8; bits1 <= FLINT_BITS; bits1 += 1)
        {
        for (bits2 = bits1; bits2 <= FLINT_BITS; bits2 += 1)
        {
            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, max_fields) + 1;
            N = words_per_exp(nfields, bits1);
            for (i = 0; i < length*nfields; i++)
                a[i] = n_randint(state, 0) & (l_shift(UWORD(1), bits1) - 1);

            /* FLINT_BITS => bits1 */
            for (i = 0; i < length; i++)
                mpoly_set_monomial(b + i*N, a + i*nfields, bits1, nfields, 0, 0);

            /* bits1 => bit2 */
            mpoly_unpack_monomials(c, bits2, b, bits1, length, nfields);

            /* bits2 => FLINT_BITS */
            mpoly_unpack_monomials(d, FLINT_BITS, c, bits2, length, nfields);

            for (i = 0; i < length*nfields; i++)
                if (a[i] != d[i])
                {
                    printf("FAIL\n");
                    flint_printf("bits1 = %wd, bits2 = %wd\n", bits1, bits2);
                    flint_abort();
                }
        }
        }
    }

    flint_free(d);
    flint_free(c);
    flint_free(b);
    flint_free(a);

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

