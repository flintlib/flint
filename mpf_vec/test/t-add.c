/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

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
#include "mpf_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("add....");
    fflush(stdout);



    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpf *a, *b, *c;
        slong len = n_randint(state, 100);

        a = _mpf_vec_init(len, 200);
        b = _mpf_vec_init(len, 200);
        c = _mpf_vec_init(len, 200);
        _mpf_vec_randtest(a, state, len, 200);
        _mpf_vec_randtest(b, state, len, 200);

        _mpf_vec_add(c, a, b, len);
        _mpf_vec_add(a, a, b, len);

        result = (_mpf_vec_equal(a, c, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            abort();
        }

        _mpf_vec_clear(a, len);
        _mpf_vec_clear(b, len);
        _mpf_vec_clear(c, len);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpf *a, *b, *c;
        slong len = n_randint(state, 100);

        a = _mpf_vec_init(len, 200);
        b = _mpf_vec_init(len, 200);
        c = _mpf_vec_init(len, 200);
        _mpf_vec_randtest(a, state, len, 200);
        _mpf_vec_randtest(b, state, len, 200);

        _mpf_vec_add(c, a, b, len);
        _mpf_vec_add(b, a, b, len);

        result = (_mpf_vec_equal(b, c, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            abort();
        }

        _mpf_vec_clear(a, len);
        _mpf_vec_clear(b, len);
        _mpf_vec_clear(c, len);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
