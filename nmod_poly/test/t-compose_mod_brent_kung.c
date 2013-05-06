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
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);
    printf("compose_mod_brent_kung....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d, e;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);
        nmod_poly_init(e, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_compose_mod_brent_kung(d, a, b, c);
        nmod_poly_compose(e, a, b);
        nmod_poly_rem(e, e, c);

        if (!nmod_poly_equal(d, e))
        {
            printf("FAIL (composition):\n");
            nmod_poly_print(a); printf("\n");
            nmod_poly_print(b); printf("\n");
            nmod_poly_print(c); printf("\n");
            nmod_poly_print(d); printf("\n");
            nmod_poly_print(e); printf("\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
        nmod_poly_clear(e);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_compose_mod_brent_kung(d, a, b, c);
        nmod_poly_compose_mod_brent_kung(a, a, b, c);

        if (!nmod_poly_equal(d, a))
        {
            printf("FAIL (aliasing a):\n");
            nmod_poly_print(a); printf("\n");
            nmod_poly_print(b); printf("\n");
            nmod_poly_print(c); printf("\n");
            nmod_poly_print(d); printf("\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_compose_mod_brent_kung(d, a, b, c);
        nmod_poly_compose_mod_brent_kung(b, a, b, c);

        if (!nmod_poly_equal(d, b))
        {
            printf("FAIL (aliasing b)\n");
            nmod_poly_print(a); printf("\n");
            nmod_poly_print(b); printf("\n");
            nmod_poly_print(c); printf("\n");
            nmod_poly_print(d); printf("\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_compose_mod_brent_kung(d, a, b, c);
        nmod_poly_compose_mod_brent_kung(c, a, b, c);

        if (!nmod_poly_equal(d, c))
        {
            printf("FAIL (aliasing c)\n");
            nmod_poly_print(a); printf("\n");
            nmod_poly_print(b); printf("\n");
            nmod_poly_print(c); printf("\n");
            nmod_poly_print(d); printf("\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
