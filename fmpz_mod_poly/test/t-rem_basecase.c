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
    flint_rand_t state;

    printf("rem_basecase....");
    fflush(stdout);

    flint_randinit(state);

    /* Check q*b + r = a, no aliasing */
    for (i = 0; i < 5000; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, q, r, t;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1);

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

        fmpz_mod_poly_divrem_basecase(q, r, a, b);
        fmpz_mod_poly_rem_basecase(t, a, b);

        result = (fmpz_mod_poly_equal(r, t));
        if (!result)
        {
            printf("FAIL (cmp with divrem):\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            printf("a = "), fmpz_mod_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_mod_poly_print(b), printf("\n\n");
            printf("q = "), fmpz_mod_poly_print(q), printf("\n\n");
            printf("r = "), fmpz_mod_poly_print(r), printf("\n\n");
            printf("t = "), fmpz_mod_poly_print(t), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_mod_poly_clear(t);
        fmpz_clear(p);
    }

    /* Alias b and r */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(r, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1);

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

        fmpz_mod_poly_rem_basecase(r, a, b);
        fmpz_mod_poly_rem_basecase(b, a, b);

        result = (fmpz_mod_poly_equal(r, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            printf("a = "), fmpz_mod_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_mod_poly_print(b), printf("\n\n");
            printf("r = "), fmpz_mod_poly_print(r), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    /* a and r */
    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(r, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1);

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

        fmpz_mod_poly_rem_basecase(r, a, b);
        fmpz_mod_poly_rem_basecase(a, a, b);

        result = (fmpz_mod_poly_equal(r, a));
        if (!result)
        {
            printf("FAIL:\n");
            printf("p = "), fmpz_print(p), printf("\n\n");
            printf("a = "), fmpz_mod_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_mod_poly_print(b), printf("\n\n");
            printf("r = "), fmpz_mod_poly_print(r), printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
