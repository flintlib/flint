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
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee

******************************************************************************/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>
#include <stdio.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("compose_mod_brent_kung_precomp_preinv....");
    fflush(stdout);

    /* no aliasing */
    for (i = 0; i < 2000; i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d, e;
        fmpz_t p;
        fmpz_mat_t B;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(cinv, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(e, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);

        fmpz_mod_poly_reverse (cinv, c, c->length);
        fmpz_mod_poly_inv_series_newton (cinv, cinv, c->length);
        fmpz_mat_init (B, n_sqrt (c->length-1)+1, c->length-1);
        fmpz_mod_poly_precompute_matrix (B, b, c, cinv);

        fmpz_mod_poly_rem(a, a, c);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(d, a, B, c, cinv);
        fmpz_mod_poly_compose(e, a, b);
        fmpz_mod_poly_rem(e, e, c);

        if (!fmpz_mod_poly_equal(d, e))
        {
            flint_printf("FAIL (composition):\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d); flint_printf("\n");
            flint_printf("e:\n"); fmpz_mod_poly_print(e); flint_printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mat_clear     (B);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(cinv);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(e);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 1000; i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d;
        fmpz_t p;
        fmpz_mat_t B;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(cinv, p);
        fmpz_mod_poly_init(d, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);

        fmpz_mod_poly_reverse (cinv, c, c->length);
        fmpz_mod_poly_inv_series_newton (cinv, cinv, c->length);
        fmpz_mat_init (B, n_sqrt (c->length-1)+1, c->length-1);
        fmpz_mod_poly_precompute_matrix (B, b, c, cinv);

        fmpz_mod_poly_rem(a, a, c);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(d, a, B, c, cinv);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(a, a, B, c, cinv);

        if (!fmpz_mod_poly_equal(d, a))
        {
            flint_printf("FAIL (aliasing a):\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d); flint_printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mat_clear     (B);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(cinv);
        fmpz_mod_poly_clear(d);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 1000; i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d;
        fmpz_t p;
        fmpz_mat_t B;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(cinv, p);
        fmpz_mod_poly_init(d, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);

        fmpz_mod_poly_reverse (cinv, c, c->length);
        fmpz_mod_poly_inv_series_newton (cinv, cinv, c->length);
        fmpz_mat_init (B, n_sqrt (c->length-1)+1, c->length-1);
        fmpz_mod_poly_precompute_matrix (B, b, c, cinv);

        fmpz_mod_poly_rem(a, a, c);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(d, a, B, c, cinv);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(c, a, B, c, cinv);

        if (!fmpz_mod_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing c)\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d); flint_printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mat_clear     (B);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(cinv);
        fmpz_mod_poly_clear(d);
    }

    /* Test aliasing of res and cinv */
    for (i = 0; i < 1000; i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d;
        fmpz_t p;
        fmpz_mat_t B;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(cinv, p);
        fmpz_mod_poly_init(d, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);

        fmpz_mod_poly_reverse (cinv, c, c->length);
        fmpz_mod_poly_inv_series_newton (cinv, cinv, c->length);
        fmpz_mat_init (B, n_sqrt (c->length-1)+1, c->length-1);
        fmpz_mod_poly_precompute_matrix (B, b, c, cinv);

        fmpz_mod_poly_rem(a, a, c);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(d, a, B, c, cinv);
        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(cinv, a, B, c, cinv);

        if (!fmpz_mod_poly_equal(d, cinv))
        {
            flint_printf("FAIL (aliasing c)\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a); flint_printf("\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b); flint_printf("\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c); flint_printf("\n");
            flint_printf("d:\n"); fmpz_mod_poly_print(d); flint_printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mat_clear     (B);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(cinv);
        fmpz_mod_poly_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
