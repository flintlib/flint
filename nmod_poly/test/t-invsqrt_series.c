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
    int i, result;
    FLINT_TEST_INIT(state);
    

    flint_printf("invsqrt_series....");
    fflush(stdout);

    /* Check 1/g^2 = h mod x^m */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t h, g, r;
        slong m;

        mp_limb_t n;
        do n = n_randtest_prime(state, 0);
        while (n == UWORD(2));

        nmod_poly_init(h, n);
        nmod_poly_init(g, n);
        nmod_poly_init(r, n);

        do nmod_poly_randtest(h, state, n_randint(state, 1000));
        while (h->length == 0);
        nmod_poly_set_coeff_ui(h, 0, UWORD(1));

        m = n_randint(state, h->length) + 1;

        nmod_poly_invsqrt_series(g, h, m);

        nmod_poly_mullow(r, g, g, m);
        nmod_poly_inv_series(r, r, m);
        nmod_poly_truncate(h, m);

        result = (nmod_poly_equal(r, h));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(h), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }
        
        nmod_poly_clear(h);
        nmod_poly_clear(g);
        nmod_poly_clear(r);
    }

    /* Check aliasing of h and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t g, h;
        slong m;

        mp_limb_t n;
        do n = n_randtest_prime(state, 0);
        while (n == UWORD(2));

        nmod_poly_init(h, n);
        nmod_poly_init(g, n);
        do nmod_poly_randtest(h, state, n_randint(state, 500));
        while (h->length == 0);
        nmod_poly_set_coeff_ui(h, 0, UWORD(1));

        m = n_randint(state, h->length) + 1;

        nmod_poly_invsqrt_series(g, h, m);
        nmod_poly_invsqrt_series(h, h, m);

        result = (nmod_poly_equal(g, h));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(h), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd, m = %wd\n", n, m);
            abort();
        }

        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
