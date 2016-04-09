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

    Copyright (C) 2016 Vincent Delecroix

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"

int
main(void)
{
    int l, result;
    mp_limb_t i, j, k, tot;
    nmod_t mod;
    mp_limb_t n;

    FLINT_TEST_INIT(state);

    flint_printf("power_sums_schoenhage....");

    /* Check that it is valid in degree 3 with integer roots, ie */
    /* for polynomials of the form (x-i)(x-j)(x-k)               */
    n = 101;                    /* TODO: a random number not divisible by small primes! */
    nmod_init(&mod, n);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {
                nmod_poly_t a, b, c;
                nmod_poly_init(a, n);
                nmod_poly_init(b, n);
                nmod_poly_init(c, n);

                nmod_poly_set_coeff_ui(a, 0, n - i * j * k);
                nmod_poly_set_coeff_ui(a, 1, i * j + i * k + j * k);
                nmod_poly_set_coeff_ui(a, 2, n - i - j - k);
                nmod_poly_set_coeff_ui(a, 3, 1);

                nmod_poly_power_sums_schoenhage(b, a, 20);

                for (l = 0; l < FLINT_MIN(20, nmod_poly_length(b)); l++)
                {
                    tot = nmod_add(nmod_pow_ui(i, l, mod),
                                   nmod_pow_ui(j, l, mod), mod);
                    tot = nmod_add(tot, nmod_pow_ui(k, l, mod), mod);

                    result = nmod_poly_get_coeff_ui(b, l) == tot;
                    if (!result)
                    {
                        flint_printf("FAIL: power sums integral root\n");
                        flint_printf("%d %d %d %d\n", i, j, k, l);
                        abort();
                    }
                }

                nmod_poly_power_sums_to_poly_schoenhage(c, b);
                result = nmod_poly_equal(a, c);
                if (!result)
                {
                    flint_printf("FAIL: power sums to poly schoenhage\n");
                    flint_printf("a = ");
                    nmod_poly_print(a), flint_printf("\n\n");
                    flint_printf("b = ");
                    nmod_poly_print(b), flint_printf("\n\n");
                    flint_printf("c = ");
                    nmod_poly_print(c), flint_printf("\n\n");
                    abort();
                }

                nmod_poly_clear(a);
                nmod_poly_clear(b);
                nmod_poly_clear(c);
            }

    /* Check that going back and forth between the power sums representation gives the identity */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        n = 5003;               /* TODO: a random number not divisible by small primes! */

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 20));
        nmod_poly_make_monic(a, a);

        nmod_poly_power_sums_schoenhage(b, a, 30);
        nmod_poly_power_sums_to_poly_schoenhage(c, b);

        result = nmod_poly_equal(a, c);
        if (!result)
        {
            flint_printf("FAIL: power sums - power sums to poly\n");
            flint_printf("a = ");
            nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = ");
            nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("c = ");
            nmod_poly_print(c), flint_printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }


    /* Check that the product of polynomials correspond to the sum of Power sums series */
    /* (and aliasing of nmod_poly_power_sums)                                           */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        n = 5003;               /* TODO: a random number not divisible by small primes! */

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 10));
        nmod_poly_randtest_not_zero(b, state, 1 + n_randint(state, 10));

        nmod_poly_mul(c, a, b);
        nmod_poly_power_sums_schoenhage(c, c, 20);

        /* NOTE: the code path is not the same if the polynomial is monic. We let only a be monic */
        nmod_poly_make_monic(a, a);
        nmod_poly_power_sums_schoenhage(a, a, 20);
        nmod_poly_power_sums_schoenhage(b, b, 20);
        nmod_poly_add(d, a, b);

        result = nmod_poly_equal(c, d);
        if (!result)
        {
            flint_printf
                ("FAIL: PowerSums(p1 p2) = PowerSums(p1) + PowerSums(p2)\n");
            flint_printf("a = ");
            nmod_poly_print(a), flint_printf("\n");
            flint_printf("b = ");
            nmod_poly_print(b), flint_printf("\n");
            flint_printf("c = ");
            nmod_poly_print(c), flint_printf("\n");
            flint_printf("d = ");
            nmod_poly_print(d), flint_printf("\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
