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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("bit_pack/bit_unpack....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 20000; i++)
    {
        fmpz_poly_t a, b;

        long length = n_randint(100, state) + 1;
        mp_bitcnt_t bits = n_randint(300, state) + 2;
        mp_ptr arr = (mp_ptr) calloc((length * bits - 1) / FLINT_BITS + 1,
                                     sizeof(mp_limb_t));
        int negate;

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        /* -1 bit to handle signs */
        fmpz_poly_randtest_not_zero(a, state, length, bits - 1);

        negate = fmpz_sgn(a->coeffs + a->length - 1);
        if (negate > 0)
            negate = 0;

        _fmpz_poly_bit_pack(arr, a->coeffs, a->length, bits, negate);
        fmpz_poly_fit_length(b, a->length);
        _fmpz_poly_bit_unpack(b->coeffs, a->length, arr, bits, negate);
        _fmpz_poly_set_length(b, a->length);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        free(arr);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    for (i = 0; i < 20000; i++)
    {
        fmpz_poly_t a, b;

        long length = n_randint(100, state) + 1;
        mp_bitcnt_t bits = n_randint(300, state) + 1;
        mp_ptr arr = (mp_ptr) calloc((length * bits - 1) / FLINT_BITS + 1,
                                     sizeof(mp_limb_t));

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        do
            fmpz_poly_randtest_unsigned(a, state, length, bits);
        while (a->length == 0);

        _fmpz_poly_bit_pack(arr, a->coeffs, a->length, bits, 0);
        fmpz_poly_fit_length(b, a->length);
        _fmpz_poly_bit_unpack_unsigned(b->coeffs, a->length, arr, bits);
        _fmpz_poly_set_length(b, a->length);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        free(arr);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
