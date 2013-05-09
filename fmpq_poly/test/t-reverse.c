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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("reverse....");
    fflush(stdout);

    flint_randinit(state);

    /* Aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        len_t n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 150);

        fmpq_poly_reverse(a, b, n);
        fmpq_poly_reverse(b, b, n);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            printf("a = "), fmpq_poly_print(a), printf("\n\n");
            printf("b = "), fmpq_poly_print(b), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Correctness (?) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        len_t j, len, n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 150);

        len = FLINT_MIN(n, b->length);
        if (len)
        {
            fmpq_poly_fit_length(a, n);
            for (j = 0; j < len; j++)
                fmpz_set(a->coeffs + (n - len) + j, b->coeffs + (len - 1 - j));
            fmpz_set(a->den, b->den);
            _fmpq_poly_set_length(a, n);
            fmpq_poly_canonicalise(a);
        }

        fmpq_poly_reverse(b, b, n);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            printf("a = "), fmpq_poly_print(a), printf("\n\n");
            printf("b = "), fmpq_poly_print(b), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
