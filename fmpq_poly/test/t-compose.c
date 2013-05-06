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
    ulong cflags = 0UL;

    printf("compose....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpq_poly_compose(f, g, h);
        fmpq_poly_compose(g, g, h);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            printf("FAIL (aliasing 1):\n");
            fmpq_poly_debug(f), printf("\n\n");
            fmpq_poly_debug(g), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpq_poly_compose(f, g, h);
        fmpq_poly_compose(h, g, h);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(h) ? 0 : 2;
        result = (fmpq_poly_equal(f, h) && !cflags);
        if (!result)
        {
            printf("FAIL (aliasing 2):\n");
            fmpq_poly_debug(f), printf("\n\n");
            fmpq_poly_debug(h), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Compare with the naive method for g(h(t)) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h, s, t, u;
        mpq_t c;
        long k;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_init(u);
        mpq_init(c);
        fmpq_poly_randtest(g, state, n_randint(state, 20), 65);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 65);
        
        fmpq_poly_zero(s);
        fmpq_poly_set_ui(t, 1UL);
        for (k = 0L; k < g->length; k++)
        {
            fmpq_poly_get_coeff_mpq(c, g, k);
            fmpq_poly_scalar_mul_mpq(u, t, c);
            fmpq_poly_add(s, s, u);
            fmpq_poly_mul(t, t, h);
        }
        
        fmpq_poly_compose(f, g, h);

        result = (fmpq_poly_equal(f, s));
        if (!result)
        {
            printf("FAIL (compare with naive):\n");
            printf("g = "), fmpq_poly_debug(g), printf("\n\n");
            printf("h = "), fmpq_poly_debug(h), printf("\n\n");
            printf("f = "), fmpq_poly_debug(f), printf("\n\n");
            printf("s = "), fmpq_poly_debug(s), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
        fmpq_poly_clear(u);
        mpq_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
