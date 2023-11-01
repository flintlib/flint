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

TEST_FUNCTION_START(mpoly_max_degrees_tight, state)
{
    slong k, i, j, length, nfields, bits1, bits2;
    slong * bases, * max, * max2, * prods;
    ulong * a, * b, * c, * t;
    ulong max_length, max_fields;

    max_length = 100;
    max_fields = FLINT_BITS/8;  /* exponents should fit in one word */

    a = flint_malloc(max_length*max_fields*sizeof(ulong));
    b = flint_malloc(max_length*max_fields*sizeof(ulong));
    c = flint_malloc(max_length*max_fields*sizeof(ulong));
    t = flint_malloc(max_length*sizeof(ulong));
    bases = (slong *) flint_malloc(max_fields*sizeof(slong));
    prods = (slong *) flint_malloc((max_fields + 1)*sizeof(slong));
    max =   (slong *) flint_malloc(max_fields*sizeof(ulong));
    max2 =  (slong *) flint_malloc(max_fields*sizeof(ulong));

    for (k = 0; k < 1000 * flint_test_multiplier(); k++)
    {
        /* do FLINT_BITS => bits1
                         => tight packing
                         => bits2 => FLINT_BITS and compare */
        for (bits1 = 8; bits1 <= FLINT_BITS; bits1 *= 2)
        for (bits2 = bits1; bits2 <= FLINT_BITS; bits2 *= 2)
        {
            length = n_randint(state, max_length) + 1;
            nfields = n_randint(state, FLINT_BITS/FLINT_MAX(bits1, bits2)) + 1;

            for (j = 0; j < nfields; j++)
                max[j] = 0;

            for (j = 0; j < nfields; j++)
                bases[j] =  n_randint(state, 200) + 1;

            prods[0] = 1;
            for (i = 0; i < nfields; i++)
                prods[i + 1] = prods[i]*bases[i];

            for (i = 0; i < nfields*length; i += nfields)
                for (j = 0; j < nfields; j++)
                {
                    a[i + j] = n_randint(state, bases[j]);
                    max[j] = FLINT_MAX(max[j], a[i + j]);
                }

            mpoly_pack_vec_ui(b, a, bits1, nfields, length);
            mpoly_pack_monomials_tight(t, b, length, bases, nfields, bits1);
            mpoly_unpack_monomials_tight(c, t, length, bases, nfields, bits2);

            mpoly_max_degrees_tight(max2, t, length, prods, nfields);

            for (j = 0; j < nfields; j++)
                if (max[j] != max2[j])
                {
                    flint_printf("FAIL\nmax_degrees_tight");
                    fflush(stdout);
                    flint_abort();
                }
        }
    }

    flint_free(max2);
    flint_free(max);
    flint_free(prods);
    flint_free(bases);
    flint_free(t);
    flint_free(c);
    flint_free(b);
    flint_free(a);

    TEST_FUNCTION_END(state);
}
