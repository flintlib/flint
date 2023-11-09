/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_bit_pack, state)
{
    int i, result;

    for (i = 0; i < 50000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        flint_bitcnt_t bits = n_randint(state, 300) + 1;
        ulong space = (300 - 1) / FLINT_BITS + 2;   /* 2 to accommodate shift */
        mp_ptr arr = (mp_ptr) flint_calloc(space, sizeof(mp_limb_t));
        flint_bitcnt_t shift = n_randint(state, FLINT_BITS);
        int negate = (int) -n_randint(state, 2);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest(a, state, bits - 1); /* need one bit for sign */

        arr[0] = n_randbits(state, shift);

        fmpz_bit_pack(arr, shift, bits, a, negate, 0);
        fmpz_bit_unpack(b, arr, shift, bits, negate, 0);

        result = (fmpz_cmp(a, b) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n");
            fmpz_print(b), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(arr);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    for (i = 0; i < 50000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        flint_bitcnt_t bits = n_randint(state, 300) + 1;
        ulong space = (300 - 1) / FLINT_BITS + 2;   /* 2 to accommodate shift */
        mp_ptr arr = (mp_ptr) flint_calloc(space, sizeof(mp_limb_t));
        flint_bitcnt_t shift = n_randint(state, FLINT_BITS);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest_unsigned(a, state, bits);

        arr[0] = n_randbits(state, shift);

        fmpz_bit_pack(arr, shift, bits, a, 0, 0);
        fmpz_bit_unpack_unsigned(b, arr, shift, bits);

        result = (fmpz_cmp(a, b) == 0) && _fmpz_is_canonical(b);

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n");
            fmpz_print(b), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(arr);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
