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

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("divides....");
    fflush(stdout);

    flint_randinit(state);

    /* Check that b divides a*b and that the quotient is a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;
        
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        fmpz_poly_mul(p, a, b);

        result = (fmpz_poly_divides(q, p, b) && fmpz_poly_equal(q, a));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            fmpz_poly_print(p), printf("\n\n");
            fmpz_poly_print(q), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    /* Check aliasing of q with a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p;
        
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        fmpz_poly_mul(p, a, b);

        result = (fmpz_poly_divides(p, p, b) && fmpz_poly_equal(p, a));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            fmpz_poly_print(p), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
    }

    /* Check aliasing of q with b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p;
        
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        fmpz_poly_mul(p, a, b);

        result = (fmpz_poly_divides(b, p, b) && fmpz_poly_equal(b, a));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            fmpz_poly_print(p), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
    }

    /* Check when not divisible */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q, g, s;
        
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_init(s);
        fmpz_poly_init(g);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        do {
           fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        } while (b->length < 2);
        fmpz_poly_mul(p, a, b);
        do {
           fmpz_poly_randtest_not_zero(s, state, b->length, 200);
           fmpz_poly_gcd(g, s, b);
        } while (g->length == b->length);
        fmpz_poly_add(p, p, s);

        result = (!fmpz_poly_divides(q, p, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            fmpz_poly_print(p), printf("\n\n");
            fmpz_poly_print(q), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
        fmpz_poly_clear(s);
        fmpz_poly_clear(g);
    }

    
    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
