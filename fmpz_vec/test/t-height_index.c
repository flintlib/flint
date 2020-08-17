/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

static slong
refimpl(const fmpz * v, slong len)
{
    slong i, max = 0;

    for (i = 1; i < len; i++)
        if (fmpz_cmpabs(v + i, v + max) > 0)
            max = i;

    return max;
}

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("height_index....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        slong len, bits, p1, p2;

        len = 1 + n_randint(state, 100);

        a = _fmpz_vec_init(len);
        bits = n_randint(state, 200);
        _fmpz_vec_randtest(a, state, len, bits);

        p1 = _fmpz_vec_height_index(a, len);
        p2 = refimpl(a, len);

        result = (p1 == p2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("bits = %wd, p1 = %wd, p2 = %wd\n", bits, p1, p2);
            abort();
        }

        _fmpz_vec_clear(a, len);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
