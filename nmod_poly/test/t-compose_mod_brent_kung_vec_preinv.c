/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013, 2014 Martin Lee

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
    
    flint_printf("compose_mod_brent_kung_vec_preinv....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, ainv, b, c;
        mp_limb_t m = n_randtest_prime(state, 0);
        slong j, k, l;
        nmod_poly_struct * pow, * res;

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(ainv, m);

        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(a, state, 1+n_randint(state, 20));
        l= n_randint(state, 20) + 1;
        k= n_randint(state, l ) + 1;

        nmod_poly_rem(b, b, a);
        nmod_poly_reverse(ainv, a, a->length);
        nmod_poly_inv_series(ainv, ainv, a->length);
        pow = (nmod_poly_struct *) flint_malloc((l + k)*sizeof(nmod_poly_struct));
        res = pow + l;

        for (j = 0; j < l; j++)
        {
            nmod_poly_init(pow + j, m);
            nmod_poly_randtest(pow + j, state, n_randint(state, 20) + 1);
            nmod_poly_rem(pow + j, pow + j, a);
        }

	for (j = 0; j < k; j++)
	    nmod_poly_init(res + j, m);

        nmod_poly_compose_mod_brent_kung_vec_preinv(res, pow, l, k, b, a, ainv);

        for (j = 0; j < k; j++)
        {
            nmod_poly_compose_mod(c, pow + j, b, a);
            if (!nmod_poly_equal(res + j, c))
            {
                flint_printf("FAIL (composition):\n");
                flint_printf("a:\n"); nmod_poly_print(a); flint_printf("\n");
                flint_printf("res:\n"); nmod_poly_print(res + j); flint_printf("\n");
                flint_printf("pow:\n"); nmod_poly_print(pow + j); flint_printf("\n");
                flint_printf("b:\n"); nmod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); nmod_poly_print(c); flint_printf("\n");
                flint_printf("j: %wd\n", j);
                abort();
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(ainv);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        for (j = 0; j < l; j++)
            nmod_poly_clear(pow + j);
        for (j = 0; j < k; j++)
            nmod_poly_clear(res + j);
        flint_free(pow);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
