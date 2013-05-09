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
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mul_karatsuba....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        fmpz_poly_mul_karatsuba(a, b, c);
        fmpz_poly_mul_karatsuba(b, b, c);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        fmpz_poly_mul_karatsuba(a, b, c);
        fmpz_poly_mul_karatsuba(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Compare with mul_classical */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        fmpz_poly_mul_karatsuba(a, b, c);
        fmpz_poly_mul_classical(d, b, c);

        result = (fmpz_poly_equal(a, d));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(d), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* Check _fmpz_poly_mul_karatsuba directly */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        len_t len1, len2;
        fmpz_poly_t a, b, out1, out2;

        len1 = n_randint(state, 100) + 1;
        len2 = n_randint(state, 100) + 1;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(out1);
        fmpz_poly_init(out2);
        fmpz_poly_randtest(a, state, len1, 200);
        fmpz_poly_randtest(b, state, len2, 200);

        fmpz_poly_mul_karatsuba(out1, a, b);
        fmpz_poly_fit_length(a, a->alloc + n_randint(state, 10));
        fmpz_poly_fit_length(b, b->alloc + n_randint(state, 10));
        a->length = a->alloc;
        b->length = b->alloc;
        fmpz_poly_fit_length(out2, a->length + b->length - 1);
        if (a->length >= b->length)
            _fmpz_poly_mul_karatsuba(out2->coeffs, a->coeffs, a->length,
                                                   b->coeffs, b->length);
        else
            _fmpz_poly_mul_karatsuba(out2->coeffs, b->coeffs, b->length,
                                                   a->coeffs, a->length);
        _fmpz_poly_set_length(out2, a->length + b->length - 1);
        _fmpz_poly_normalise(out2);

        result = (fmpz_poly_equal(out1, out2));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(out1), printf("\n\n");
            fmpz_poly_print(out2), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(out1);
        fmpz_poly_clear(out2);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
