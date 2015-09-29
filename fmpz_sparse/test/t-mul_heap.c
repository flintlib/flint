/*============================================================================

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

    Copyright (C) 509 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("mul_heap....");
    fflush(stdout);

    
    /* Check aliasing of a and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a, b, c;
        fmpz_t d, e;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest(d, state, 50);
        fmpz_randtest(e, state, 50);

        fmpz_sparse_init(a);
        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_randtest(b, state, n_randint(state, 50), d, 50);
        fmpz_sparse_randtest(c, state, n_randint(state, 50), e, 50);

        fmpz_sparse_mul_heaps(a, b, c);
        fmpz_sparse_mul_heaps(b, b, c);

        result = (fmpz_sparse_equal(a, b));
        if (!result)
        {
          flint_printf("FAIL PHASE 1:\n");
          fmpz_sparse_print(a), flint_printf("\n\n");
          fmpz_sparse_print(b), flint_printf("\n\n");
          abort();
        }

        fmpz_sparse_clear(a);
        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a, b, c;
        fmpz_t d, e;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest(d, state, 50);
        fmpz_randtest(e, state, 50);

        fmpz_sparse_init(a);
        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_randtest(b, state, n_randint(state, 50), d, 50);
        fmpz_sparse_randtest(c, state, n_randint(state, 50), e, 50);

        fmpz_sparse_mul_heaps(a, b, c);
        fmpz_sparse_mul_heaps(c, b, c);

        result = (fmpz_sparse_equal(a, c));
        if (!result)
        {
          flint_printf("FAIL PHASE 2:\n");
          fmpz_sparse_print(a), flint_printf("\n\n");
          fmpz_sparse_print(c), flint_printf("\n\n");
          abort();
        }

        fmpz_sparse_clear(a);
        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a1, a2, b, c, d;
        fmpz_t e, f, g;

        fmpz_init(e);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_randtest(e, state, 50);
        fmpz_randtest(f, state, 50);
        fmpz_randtest(g, state, 50);

        fmpz_sparse_init(a1);
        fmpz_sparse_init(a2);
        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_init(d);
        fmpz_sparse_randtest(b, state, n_randint(state, 50), e, 50);
        fmpz_sparse_randtest(c, state, n_randint(state, 50), f, 50);
        fmpz_sparse_randtest(d, state, n_randint(state, 50), g, 50);

        fmpz_sparse_mul_heaps(a1, b, c);
        fmpz_sparse_mul_heaps(a2, b, d);
        fmpz_sparse_new_add(a1, a1, a2);

        fmpz_sparse_new_add(c, c, d);
        fmpz_sparse_mul_heaps(a2, b, c);

        result = (fmpz_sparse_equal(a1, a2));
        if (!result)
        {
          flint_printf("FAIL PHASE 3:\n");
          fmpz_sparse_print(a1), flint_printf("\n\n");
          fmpz_sparse_print(a2), flint_printf("\n\n");
          abort();
        }

        fmpz_sparse_clear(a1);
        fmpz_sparse_clear(a2);
        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        fmpz_sparse_clear(d);
        fmpz_clear(e);
        fmpz_clear(f);
        fmpz_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
