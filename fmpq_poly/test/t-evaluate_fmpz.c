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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("evaluate_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check that (f+g)(a) = f(a) + g(a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        fmpq_poly_t f, g, h;
        fmpq_t x, y;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(a);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);
        fmpz_randtest(a, state, n_randint(state, 100));

        fmpq_poly_evaluate_fmpz(x, f, a);
        fmpq_poly_evaluate_fmpz(y, g, a);
        fmpq_add(x, x, y);
        fmpq_poly_add(h, f, g);
        fmpq_poly_evaluate_fmpz(y, h, a);

        result = (fmpq_equal(x, y));
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpq_poly_debug(f), printf("\n");
            printf("g = "), fmpq_poly_debug(g), printf("\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("f(a) + g(a) = "), fmpq_print(x), printf("\n\n");
            printf("(f + g)(a)  = "), fmpq_print(y), printf("\n\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(a);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check that (f*g)(a) = f(a) * g(a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        fmpq_poly_t f, g;
        fmpq_t x, y;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(a);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);
        fmpz_randtest(a, state, n_randint(state, 100));

        fmpq_poly_evaluate_fmpz(x, f, a);
        fmpq_poly_evaluate_fmpz(y, g, a);
        fmpq_mul(x, x, y);
        fmpq_poly_mul(f, f, g);
        fmpq_poly_evaluate_fmpz(y, f, a);

        result = (fmpq_equal(x, y));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(a);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
