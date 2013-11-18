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

    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("derivative....");
    fflush(stdout);

    

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_set(b, a);

        fmpz_mod_poly_derivative(c, b);
        fmpz_mod_poly_derivative(b, b);

        result = (fmpz_mod_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (alias):\n");
            fmpz_mod_poly_print(a), flint_printf("\n\n");
            fmpz_mod_poly_print(b), flint_printf("\n\n");
            fmpz_mod_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_clear(p);
    }

    /* Check constants have derivative zero */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 2));

        fmpz_mod_poly_derivative(b, a);

        result = (b->length == 0);
        if (!result)
        {
            flint_printf("FAIL (da == 0):\n");
            fmpz_mod_poly_print(a), flint_printf("\n\n");
            fmpz_mod_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_clear(p);
    }

    /* Check (f g)' = f' g + f g' */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d, lhs, rhs;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(lhs, p);
        fmpz_mod_poly_init(rhs, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));

        fmpz_mod_poly_mul(lhs, a, b);
        fmpz_mod_poly_derivative(lhs, lhs);
        fmpz_mod_poly_derivative(c, a);
        fmpz_mod_poly_derivative(d, b);
        fmpz_mod_poly_mul(c, c, b);
        fmpz_mod_poly_mul(d, a, d);
        fmpz_mod_poly_add(rhs, c, d);

        result = fmpz_mod_poly_equal(lhs, rhs);
        if (!result)
        {
            flint_printf("FAIL (Leibniz):\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("(ab)' = "), fmpz_mod_poly_print(lhs), flint_printf("\n\n");
            flint_printf("a'b + ab' = "), fmpz_mod_poly_print(rhs), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(lhs);
        fmpz_mod_poly_clear(rhs);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

