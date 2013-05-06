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
#include "fmpz.h"
#include "fmpq.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("cmp....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        mpq_t X, Y;
        int c1, c2;

        fmpq_init(x);
        fmpq_init(y);
        mpq_init(X);
        mpq_init(Y);

        fmpq_randtest(x, state, 200);
        fmpq_randtest(y, state, 200);

        fmpq_get_mpq(X, x);
        fmpq_get_mpq(Y, y);

        c1 = fmpq_cmp(x, y);
        c2 = mpq_cmp(X, Y);

        if (c1 < 0) c1 = -1;
        if (c1 > 0) c1 = 1;

        if (c2 < 0) c2 = -1;
        if (c2 > 0) c2 = 1;

        if (c1 != c2)
        {
            printf("FAIL\n");
            printf("x = ");
            fmpq_print(x);
            printf("\ny = ");
            fmpq_print(y);
            printf("\ncmp(x,y) = %d, cmp(X,Y) = %d\n", c1, c2);
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

