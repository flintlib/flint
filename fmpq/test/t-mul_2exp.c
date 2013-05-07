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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("mul_2exp....");
    fflush(stdout);

    /* x = y * 2^exp */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;
        mp_bitcnt_t c;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_numref(y), fmpq_numref(y), n_randint(state, 200));
        else if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), n_randint(state, 200));
        fmpq_canonicalise(y);

        c = n_randint(state, 200);
        fmpq_mul_2exp(x, y, c);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical!\n");
            printf("x = ");
            fmpq_print(x);
            printf("\ny = ");
            fmpq_print(y);
            printf("\nc = %lu\n", c);
            abort();
        }

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);
        mpq_mul_2exp(X, Y, c);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            printf("FAIL: fmpq_mul_2exp(x,y,c) != mpq_mul_2exp(X,Y,c)\n");
            printf("x = ");
            fmpq_print(x);
            printf("\ny = ");
            fmpq_print(y);
            printf("\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        mpq_clear(X);
        mpq_clear(Y);
    }

    /* y = y * 2^exp */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mp_bitcnt_t c;

        fmpq_init(x);
        fmpq_init(y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_numref(y), fmpq_numref(y), n_randint(state, 200));
        else if (n_randint(state, 5) == 0)
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), n_randint(state, 200));
        fmpq_canonicalise(y);

        c = n_randint(state, 200);
        fmpq_mul_2exp(x, y, c);
        fmpq_mul_2exp(y, y, c);

        if (!fmpq_is_canonical(y))
        {
            printf("FAIL: result not canonical!\n");
            printf("x = ");
            fmpq_print(x);
            printf("\ny = ");
            fmpq_print(y);
            printf("\nc = %lu\n", c);
            abort();
        }

        if (!fmpq_equal(x, y))
        {
            printf("FAIL: fmpq_mul_2exp(x,y,c) != fmpq_mul_2exp(y,y,c)\n");
            printf("x = ");
            fmpq_print(x);
            printf("\ny = ");
            fmpq_print(y);
            printf("\nc = %lu\n", c);
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
