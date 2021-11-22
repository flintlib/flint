/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>
#include <stdio.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("compose_mod_brent_kung_preinv....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, cinv, d, e;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(cinv, m);
        nmod_poly_init(d, m);
        nmod_poly_init(e, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);
        nmod_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv);
        nmod_poly_compose(e, a, b);
        nmod_poly_rem(e, e, c);

        if (!nmod_poly_equal(d, e))
        {
            flint_printf("FAIL (composition):\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(cinv); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            nmod_poly_print(e); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);
        nmod_poly_clear(e);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, cinv, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(cinv, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);
        nmod_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv);
        nmod_poly_compose_mod_brent_kung_preinv(a, a, b, c, cinv);

        if (!nmod_poly_equal(d, a))
        {
            flint_printf("FAIL (aliasing a):\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(cinv); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, cinv, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(cinv, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);
        nmod_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv);
        nmod_poly_compose_mod_brent_kung_preinv(b, a, b, c, cinv);

        if (!nmod_poly_equal(d, b))
        {
            flint_printf("FAIL (aliasing b)\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(cinv); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, cinv, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(cinv, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);
        nmod_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv);
        nmod_poly_compose_mod_brent_kung_preinv(c, a, b, c, cinv);

        if (!nmod_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing c)\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(cinv); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);
    }

    /* Test aliasing of res and cinv */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, cinv, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(cinv, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_rem(a, a, c);
        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);
        nmod_poly_compose_mod_brent_kung_preinv(d, a, b, c, cinv);
        nmod_poly_compose_mod_brent_kung_preinv(cinv, a, b, c, cinv);

        if (!nmod_poly_equal(d, cinv))
        {
            flint_printf("FAIL (aliasing cinv)\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(cinv); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
