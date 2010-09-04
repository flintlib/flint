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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    fmpz_randstate_t state;

    printf("pow....");
    fflush(stdout);

    fmpq_poly_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 2000; i++)
    {
        fmpq_poly_t a, b;
        ulong exp;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(b, state, n_randint(10), 100);

        exp = (ulong) n_randtest() % 20UL;

        fmpq_poly_pow(a, b, exp);
        fmpq_poly_pow(b, b, exp);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("exp = %lu\n", exp);
            printf("a = "), fmpq_poly_print(a), printf("\n\n");
            printf("b = "), fmpq_poly_print(b), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Compare with repeated multiplications by the base */
    for (i = 0; i < 2000; i++)
    {
        fmpq_poly_t a, b, c;
        ulong exp;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(b, state, n_randint(10), 100);

        exp = (ulong) n_randtest() % 20UL;

        fmpq_poly_pow(a, b, exp);

        if (exp == 0UL && b->length > 0)
        {
            fmpq_poly_fit_length(c, 1);
            fmpz_set_ui(c->coeffs, 1UL);
            _fmpq_poly_set_length(c, 1);
        }
        else
        {
            ulong j;
            fmpq_poly_set(c, b);

            for (j = 1; j < exp; j++)
                fmpq_poly_mul(c, c, b);
        }

        result = (fmpq_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            printf("exp = %lu\n", exp);
            printf("a = "), fmpq_poly_print(a), printf("\n\n");
            printf("c = "), fmpq_poly_print(c), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    fmpq_poly_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
