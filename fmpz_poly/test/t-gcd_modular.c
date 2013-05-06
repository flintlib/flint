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
#include "fmpz_poly.h"
#include "ulong_extras.h"

/*
    Tests whether the polynomial is suitably normalised for the 
    result of a GCD operation, that is, whether it's leading 
    coefficient is non-negative.
 */
static 
int _t_gcd_is_canonical(const fmpz_poly_t poly)
{
    return fmpz_poly_is_zero(poly) || (fmpz_sgn(fmpz_poly_lead(poly)) > 0);
}

int
main(void)
{
    int i, result;
    flint_rand_t state;
    
    printf("gcd_modular....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd_modular(a, b, c);
        fmpz_poly_gcd_modular(b, b, c);

        result = (fmpz_poly_equal(a, b) && _t_gcd_is_canonical(a));
        if (!result)
        {
            printf("FAIL (aliasing a and b):\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd_modular(a, b, c);
        fmpz_poly_gcd_modular(c, b, c);

        result = (fmpz_poly_equal(a, c) && _t_gcd_is_canonical(a));
        if (!result)
        {
            printf("FAIL (aliasing a and c):\n");
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("c = "), fmpz_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check that a divides GCD(af, ag) */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 100) + 1, 40);
        fmpz_poly_randtest(f, state, n_randint(state, 100), 40);
        fmpz_poly_randtest(g, state, n_randint(state, 100), 40);

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        fmpz_poly_gcd_modular(d, f, g);

        fmpz_poly_divrem_divconquer(q, r, d, a);

        result = fmpz_poly_is_zero(r) && _t_gcd_is_canonical(d);
        if (!result)
        {
           printf("FAIL (check a | gcd(af, ag)):\n");
           printf("f = "), fmpz_poly_print(f), printf("\n");
           printf("g = "), fmpz_poly_print(g), printf("\n");
           printf("a = "), fmpz_poly_print(a), printf("\n");
           printf("d = "), fmpz_poly_print(d), printf("\n");
           abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Check that a == GCD(af, ag) when GCD(f, g) = 1 */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 100) + 1, 200);
        do {
           fmpz_poly_randtest(f, state, n_randint(state, 100), 200);
           fmpz_poly_randtest(g, state, n_randint(state, 100), 200);
           fmpz_poly_gcd_heuristic(d, f, g);
        } while (!(d->length == 1 && fmpz_is_one(d->coeffs)));

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        fmpz_poly_gcd_modular(d, f, g);

        if (!_t_gcd_is_canonical(a)) fmpz_poly_neg(a, a);

        result = fmpz_poly_equal(d, a) && _t_gcd_is_canonical(d);
        if (!result)
        {
           printf("FAIL (check a == gcd(af, ag) when gcd(f, g) = 1):\n");
           printf("f = "), fmpz_poly_print(f), printf("\n");
           printf("g = "), fmpz_poly_print(g), printf("\n");
           printf("a = "), fmpz_poly_print(a), printf("\n");
           printf("d = "), fmpz_poly_print(d), printf("\n");
           abort();
        } 

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Sebastian's test case */
    {
        fmpz_poly_t a, b, d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(d);

        fmpz_poly_set_coeff_ui(b, 2, 1);
        fmpz_poly_set_coeff_si(a, 0, -32);
        fmpz_poly_set_coeff_si(a, 1, 24);

        fmpz_poly_gcd_modular(d, a, b);

        result = (d->length == 1 && fmpz_is_one(d->coeffs));
        if (!result)
        {
            printf("FAIL (check 1 == gcd(x^2, 24*x - 32):\n");
            fmpz_poly_print(d); printf("\n"); 
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(d);
    }

    /* another test case */
    {
        fmpz_poly_t a, b, d, e;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(d);
        fmpz_poly_init(e);

        fmpz_poly_set_str(a, "12  0 0 0 0 0 0 0 0 0 8582594367 -9297159048333985579007 33822867456");
        fmpz_poly_set_str(b, "8  0 0 -258272396248218664896 0 -2762 -549690802047 -3771028 8796059467776");
        fmpz_poly_set_str(e, "3  0 0 1");

        fmpz_poly_gcd_modular(d, a, b);

        result = fmpz_poly_equal(d, e);
        if (!result)
        {
            printf("FAIL (check special #2):\n");
            fmpz_poly_print(d); printf("\n"); 
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(d);
        fmpz_poly_clear(e);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
