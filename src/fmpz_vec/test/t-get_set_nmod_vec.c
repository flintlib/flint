/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_get_set_nmod_vec, state)
{
    int i, result;

    /* Check conversion to and from nmod_vec */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        mp_ptr c;
        nmod_t mod;
        slong i;
        mp_limb_t t;

        slong len = n_randint(state, 100);
        mp_limb_t n = n_randtest_not_zero(state);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _nmod_vec_init(len);

        nmod_init(&mod, n);

        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_get_nmod_vec(c, a, len, mod);
        _fmpz_vec_set_nmod_vec(b, c, len, mod);

        for (i = 0; i < len; i++)
        {
            fmpz_mod_ui(a + i, a + i, n);
            t = fmpz_get_ui(a + i);
            if (t > n / 2)
                fmpz_sub_ui(a + i, a + i, n);
        }

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _nmod_vec_clear(c);
    }

    TEST_FUNCTION_END(state);
}
