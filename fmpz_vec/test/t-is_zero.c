/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("is_zero....");
    fflush(stdout);

    

    /* Check zero vector */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        slong len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_zero(a, len);

        result = (_fmpz_vec_is_zero(a, len));
        if (!result)
        {
            flint_printf("FAIL1:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
    }

    /* Check non-zero vector */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        slong len = n_randint(state, 100) + 1;

        a = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        fmpz_set_ui(a + (len - 1), UWORD(1));

        result = (!_fmpz_vec_is_zero(a, len));
        if (!result)
        {
            flint_printf("FAIL2:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
