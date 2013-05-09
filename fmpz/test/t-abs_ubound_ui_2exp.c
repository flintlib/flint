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

    Copyright (C) 2012 Fredrik Johansson

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
    len_t iter;
    int result;
    flint_rand_t state;

    printf("abs_ubound_ui_2exp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_t x, y;
        len_t bits, yexp;
        len_t exp;
        mp_limb_t man;

        fmpz_init(x);
        fmpz_init(y);

        fmpz_randtest_not_zero(x, state, 1 + n_randint(state, 400));

        bits = 1 + n_randint(state, FLINT_BITS - 1);

        /* compute an exactly rounded mantissa */
        fmpz_abs(y, x);

        if (fmpz_is_zero(y))
        {
            yexp = 0;
        }
        else
        {
            yexp = fmpz_bits(y) - bits;

            if (yexp >= 0)
            {
                fmpz_cdiv_q_2exp(y, y, yexp);
                if (fmpz_bits(y) == bits + 1)
                {
                    fmpz_tdiv_q_2exp(y, y, 1);
                    yexp--;
                }
            }
            else
            {
                fmpz_mul_2exp(y, y, -yexp);
            }
        }

        man = fmpz_abs_ubound_ui_2exp(&exp, x, bits);

        if (FLINT_BIT_COUNT(man) != bits)
        {
            printf("wrong number of bits!\n");
            printf("bits = %ld, count = %u\n\n", bits, FLINT_BIT_COUNT(man));
            printf("x = "); fmpz_print(x); printf("\n\n");
            printf("bits(x) = %ld\n\n", fmpz_bits(x));
            printf("y = "); fmpz_print(y); printf("\n\n");
            printf("yexp = %ld\n\n", yexp);
            printf("man = %lu, exp = %ld\n", man, exp);
            abort();
        }

        /* ok if equal */
        result = (fmpz_cmp_ui(y, man) == 0);

        /* ok if mantissa is 1 larger */
        if (!result)
        {
            result = ((exp == yexp) && (fmpz_cmp_ui(y, man - 1) == 0));
        }

        /* ok if the exact mantissa is 2^r-1 and overflow to 2^r happened */
        if (!result)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_set_ui(t, man);
            fmpz_mul_ui(t, t, 2);
            fmpz_sub_ui(t, t, 1);
            result = (exp == yexp + 1) && fmpz_equal(t, y);
            fmpz_clear(t);
        }

        if (!result)
        {
            printf("different from exact ceiling division\n");
            printf("bits = %ld\n\n", bits);
            printf("x = "); fmpz_print(x); printf("\n\n");
            printf("bits(x) = %ld\n\n", fmpz_bits(x));
            printf("y = "); fmpz_print(y); printf(", yexp = %ld\n\n", yexp);
            printf("man = %lu, exp = %ld\n", man, exp);
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
