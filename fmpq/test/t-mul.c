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

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("mul....");
    fflush(stdout);

    /* x = y * z */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, z;
        mpq_t X, Y, Z;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        mpq_init(X);
        mpq_init(Y);
        mpq_init(Z);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);
        fmpq_randtest(z, state, 200);

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);
        fmpq_get_mpq(Z, z);

        fmpq_mul(x, y, z);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical!\n");
            abort();
        }

        mpq_mul(X, Y, Z);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            printf("FAIL: fmpq_mul(x,y,z) != mpq_mul(X,Y,Z)\n");
            printf("x = ");
            fmpq_print(x);
            printf("\ny = ");
            fmpq_print(y);
            printf("\nz = ");
            fmpq_print(z);
            printf("\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);

        mpq_clear(X);
        mpq_clear(Y);
        mpq_clear(Z);
    }

    /* x = x * y */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);

        fmpq_mul(x, x, y);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical!\n");
            abort();
        }

        mpq_mul(X, X, Y);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            printf("FAIL: fmpq_mul(x,x,y) != mpq_mul(X,X,Y)\n");
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

    /* x = y * x */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);

        fmpq_mul(x, y, x);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical!\n");
            abort();
        }

        mpq_mul(X, Y, X);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            printf("FAIL: fmpq_mul(x,y,x) != mpq_mul(X,Y,X)\n");
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

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
