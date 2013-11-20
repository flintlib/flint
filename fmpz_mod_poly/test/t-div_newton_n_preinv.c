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

    Copyright (C) 2011 Sebastian Pancratz
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

    flint_printf("div_newton_n_preinv....");
    fflush(stdout);

    

    /* Compare to div_basecase */
    for (i = 0; i < 5000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, q2, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(q2, p);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);
        fmpz_mod_poly_div_basecase(q, a, b);
        fmpz_mod_poly_div_newton_n_preinv(q2, a, b, binv);

        result = (fmpz_mod_poly_equal(q, q2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q), flint_printf("\n\n");
            flint_printf("q2 = "), fmpz_mod_poly_print(q2), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(q2);
        fmpz_clear(p);
    }

    /* Alias a and q */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);
        fmpz_mod_poly_div_newton_n_preinv(q, a, b, binv);
        fmpz_mod_poly_div_newton_n_preinv(a, a, b, binv);

        result = (fmpz_mod_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_clear(p);
    }

    /* Alias b and q */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);
        fmpz_mod_poly_div_newton_n_preinv(q, a, b, binv);
        fmpz_mod_poly_div_newton_n_preinv(b, a, b, binv);

        result = (fmpz_mod_poly_equal(q, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_clear(p);
    }


    /* Alias binv and q */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, binv;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);
        fmpz_mod_poly_div_newton_n_preinv(q, a, b, binv);
        fmpz_mod_poly_div_newton_n_preinv(binv, a, b, binv);

        result = (fmpz_mod_poly_equal(q, binv));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("binv = "), fmpz_mod_poly_print(binv), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
