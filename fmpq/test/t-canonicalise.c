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

    printf("canonicalise....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        fmpz_t mult;

        fmpq_init(x);
        fmpq_randtest(x, state, 200);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: expected fmpq_randtest output to be canonical\n");
            fmpq_print(x);
            printf("\n");
            abort();
        }

        fmpz_init(mult);
        fmpz_randtest_not_zero(mult, state, 200);
        fmpz_add_ui(mult, mult, 1UL);

        fmpz_mul(&x->num, &x->num, mult);
        fmpz_mul(&x->den, &x->den, mult);

        if (fmpq_is_canonical(x))
        {
            printf("FAIL: expected fmpq_is_canonical to detect common factor\n");
            fmpq_print(x);
            printf("\n");
            abort();
        }

        fmpq_canonicalise(x);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical after calling fmpq_canonicalise\n");
            fmpq_print(x);
            printf("\n");
            abort();
        }

        fmpz_neg(&x->den, &x->den);

        if (fmpq_is_canonical(x))
        {
            printf("FAIL: negative denominator reported as being canonical\n");
            fmpq_print(x);
            printf("\n");
            abort();
        }

        fmpq_canonicalise(x);

        if (!fmpq_is_canonical(x))
        {
            printf("FAIL: result not canonical after calling fmpq_canonicalise\n");
            fmpq_print(x);
            printf("\n");
            abort();
        }

        fmpz_clear(mult);
        fmpq_clear(x);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
