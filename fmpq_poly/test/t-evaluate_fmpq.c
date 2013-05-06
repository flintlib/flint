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

    printf("evaluate_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y;
        fmpq_poly_t f;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_poly_init(f);
        fmpq_poly_randtest(f, state, n_randint(state, 80), 100);
        fmpz_randtest(a, state, 80);
        fmpz_randtest_not_zero(b, state, 80);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), a);
        fmpq_canonicalise(x);

        fmpq_poly_evaluate_fmpq(y, f, x);
        fmpq_poly_evaluate_fmpq(x, f, x);

        result = (fmpq_equal(x, y));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n\n");
            fmpz_print(b), printf("\n\n");
            fmpq_poly_debug(f), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_poly_clear(f);
    }

    /* Check that (f+g)(a) = f(a) + g(a) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y, z;
        fmpq_poly_t f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 80), 100);
        fmpq_poly_randtest(g, state, n_randint(state, 80), 100);
        fmpz_randtest(a, state, 80);
        fmpz_randtest_not_zero(b, state, 80);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), a);
        fmpq_canonicalise(x);

        fmpq_poly_evaluate_fmpq(y, f, x);
        fmpq_poly_evaluate_fmpq(z, g, x);
        fmpq_add(y, y, z);
        fmpq_poly_add(f, f, g);
        fmpq_poly_evaluate_fmpq(z, f, x);

        result = (fmpq_equal(y, z));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n\n");
            fmpz_print(b), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* Check that (f*g)(a) = f(a) * g(a) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y, z;
        fmpq_poly_t f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 50), 80);
        fmpq_poly_randtest(g, state, n_randint(state, 50), 80);
        fmpz_randtest(a, state, 80);
        fmpz_randtest_not_zero(b, state, 80);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), a);
        fmpq_canonicalise(x);

        fmpq_poly_evaluate_fmpq(y, f, x);
        fmpq_poly_evaluate_fmpq(z, g, x);
        fmpq_mul(y, y, z);
        fmpq_poly_mul(f, f, g);
        fmpq_poly_evaluate_fmpq(z, f, x);

        result = (fmpq_equal(y, z));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n\n");
            fmpz_print(b), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
