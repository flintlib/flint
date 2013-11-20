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

    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart
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
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("mulmod_preinv....");
    fflush(stdout);

    

    /* Check aliasing of res and a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, res, f, finv;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(finv, p);
        fmpz_mod_poly_init(res, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);
        if (a->length >= f->length)
          fmpz_mod_poly_rem (a, a, f);
        if (b->length >= f->length)
          fmpz_mod_poly_rem (b, b, f);

        fmpz_mod_poly_reverse (finv, f, f->length);
        fmpz_mod_poly_inv_series_newton (finv, finv, f->length);

        fmpz_mod_poly_mulmod_preinv(res, a, b, f, finv);
        fmpz_mod_poly_mulmod_preinv(a, a, b, f, finv);

        result = (fmpz_mod_poly_equal(res, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(finv);
        fmpz_mod_poly_clear(res);
        fmpz_clear(p);
    }

    /* Check aliasing of res and b */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, finv, res;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(finv, p);
        fmpz_mod_poly_init(res, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);
        if (a->length >= f->length)
          fmpz_mod_poly_rem (a, a, f);
        if (b->length >= f->length)
          fmpz_mod_poly_rem (b, b, f);

        fmpz_mod_poly_reverse (finv, f, f->length);
        fmpz_mod_poly_inv_series_newton (finv, finv, f->length);

        fmpz_mod_poly_mulmod_preinv(res, a, b, f, finv);
        fmpz_mod_poly_mulmod_preinv(b, a, b, f, finv);

        result = (fmpz_mod_poly_equal(res, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(finv);
        fmpz_mod_poly_clear(res);
        fmpz_clear(p);
    }

    /* Check aliasing of res and f */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, finv, res;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(finv, p);
        fmpz_mod_poly_init(res, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);
        if (a->length >= f->length)
          fmpz_mod_poly_rem (a, a, f);
        if (b->length >= f->length)
          fmpz_mod_poly_rem (b, b, f);

        fmpz_mod_poly_reverse (finv, f, f->length);
        fmpz_mod_poly_inv_series_newton (finv, finv, f->length);

        fmpz_mod_poly_mulmod_preinv(res, a, b, f, finv);
        fmpz_mod_poly_mulmod_preinv(f, a, b, f, finv);

        result = (fmpz_mod_poly_equal(res, f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(finv);
        fmpz_mod_poly_clear(res);
        fmpz_clear(p);
    }

    /* Check aliasing of res and finv */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, finv, res;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(finv, p);
        fmpz_mod_poly_init(res, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);
        if (a->length >= f->length)
          fmpz_mod_poly_rem (a, a, f);
        if (b->length >= f->length)
          fmpz_mod_poly_rem (b, b, f);

        fmpz_mod_poly_reverse (finv, f, f->length);
        fmpz_mod_poly_inv_series_newton (finv, finv, f->length);

        fmpz_mod_poly_mulmod_preinv(res, a, b, f, finv);
        fmpz_mod_poly_mulmod_preinv(finv, a, b, f, finv);

        result = (fmpz_mod_poly_equal(res, finv));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("finv:\n"); fmpz_mod_poly_print(finv), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(finv);
        fmpz_mod_poly_clear(res);
        fmpz_clear(p);
    }

    /* No aliasing */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, res1, res2, f, finv;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(finv, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1);
        if (a->length >= f->length)
          fmpz_mod_poly_rem (a, a, f);
        if (b->length >= f->length)
          fmpz_mod_poly_rem (b, b, f);

        fmpz_mod_poly_reverse (finv, f, f->length);
        fmpz_mod_poly_inv_series_newton (finv, finv, f->length);

        fmpz_mod_poly_init(res1, p);
        fmpz_mod_poly_init(res2, p);
        fmpz_mod_poly_mulmod(res1, a, b, f);
        fmpz_mod_poly_mulmod_preinv(res2, a, b, f, finv);

        result = (fmpz_mod_poly_equal(res1, res2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); fmpz_mod_poly_print(res1), flint_printf("\n\n");
            flint_printf("res2:\n"); fmpz_mod_poly_print(res2), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(finv);
        fmpz_mod_poly_clear(res1);
        fmpz_mod_poly_clear(res2);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
