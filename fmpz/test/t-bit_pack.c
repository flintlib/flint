/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("bit_pack/bit_unpack....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 50000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mp_bitcnt_t bits = n_randint(state, 300) + 1;
        ulong space = (300 - 1) / FLINT_BITS + 2;   /* 2 to accomodate shift */
        mp_ptr arr = (mp_ptr) flint_calloc(space, sizeof(mp_limb_t));
        mp_bitcnt_t shift = n_randint(state, FLINT_BITS);
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
            printf("FAIL:\n");
            fmpz_print(a), printf("\n");
            fmpz_print(b), printf("\n");
            abort();
        }

        flint_free(arr);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    for (i = 0; i < 50000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mp_bitcnt_t bits = n_randint(state, 300) + 1;
        ulong space = (300 - 1) / FLINT_BITS + 2;   /* 2 to accomodate shift */
        mp_ptr arr = (mp_ptr) flint_calloc(space, sizeof(mp_limb_t));
        mp_bitcnt_t shift = n_randint(state, FLINT_BITS);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest_unsigned(a, state, bits);

        arr[0] = n_randbits(state, shift);

        fmpz_bit_pack(arr, shift, bits, a, 0, 0);
        fmpz_bit_unpack_unsigned(b, arr, shift, bits);

        result = (fmpz_cmp(a, b) == 0);

        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n");
            fmpz_print(b), printf("\n");
            abort();
        }

        flint_free(arr);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
