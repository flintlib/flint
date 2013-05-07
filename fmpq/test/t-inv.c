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

    printf("inv....");
    fflush(stdout);

    /* x = y * z */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, z;
        mpq_t X, Y, Z, YY, ZZ;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        mpq_init(X);
        mpq_init(Y);
        mpq_init(Z);
        mpq_init(YY);
        mpq_init(ZZ);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);
        do { fmpq_randtest(z, state, 200); } while (fmpq_is_zero(z));

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);
        fmpq_get_mpq(Z, z);

        fmpq_inv(y, z);
        fmpq_inv(x, y);

        if (!fmpq_equal(x, z))
        {
            printf("FAIL: applying inv twice did not give back the input!\n");
            abort();
        }

        if (!fmpq_is_canonical(y) || !fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical!\n");
            abort();
        }

        mpq_inv(Y, Z);
        mpq_inv(X, Z);

        fmpq_get_mpq(YY, y);
        fmpq_get_mpq(ZZ, z);

        if (!mpq_equal(Y, YY) || !mpq_equal(Z, ZZ))
        {
            printf("FAIL: fmpq_inv != mpq_inv\n");
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
        mpq_clear(YY);
        mpq_clear(ZZ);
    }

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        mpq_t X, Y;

        fmpq_init(x);
        mpq_init(X);
        mpq_init(Y);

        do { fmpq_randtest(x, state, 200); } while (fmpq_is_zero(x));

        fmpq_get_mpq(X, x);

        fmpq_inv(x, x);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical!\n");
            abort();
        }

        mpq_inv(X, X);
        fmpq_get_mpq(Y, x);

        if (!mpq_equal(X, Y))
        {
            printf("FAIL: fmpq_mul(x,x,y) != mpq_mul(X,X,Y)\n");
            printf("x = ");
            fmpq_print(x);
            printf("\n");
            abort();
        }

        fmpq_clear(x);

        mpq_clear(X);
        mpq_clear(Y);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
