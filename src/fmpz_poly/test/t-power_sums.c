/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_power_sums, state)
{
    int i, j, k, l, result;
    fmpz_t il, jl, kl, tot;
    slong n;
    fmpz_poly_t a, b, c, d, e, f;

    /* Check that it is valid in degree 3 with integer roots, ie */
    /* for polynomials of the form (x-i)(x-j)(x-k)               */
    for (i = -5; i < 5; i++)
        for (j = i; j < 5; j++)
            for (k = j; k < 5; k++)
            {
                fmpz_poly_init(a);
                fmpz_poly_init(b);
                fmpz_poly_init(c);

                fmpz_poly_set_coeff_si(a, 0, -i * j * k);
                fmpz_poly_set_coeff_si(a, 1, i * j + i * k + j * k);
                fmpz_poly_set_coeff_si(a, 2, -i - j - k);
                fmpz_poly_set_coeff_si(a, 3, 1);

                fmpz_poly_power_sums(b, a, 20);
                fmpz_poly_power_sums_naive(c, a, 20);
                result = fmpz_poly_equal(b, c);
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("%d %d %d\n\n", i, j, k);
                    fmpz_poly_print(b), flint_printf("\n\n");
                    fmpz_poly_print(c), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_init(il);
                fmpz_init(jl);
                fmpz_init(kl);
                fmpz_set_ui(il, 1);
                fmpz_set_ui(jl, 1);
                fmpz_set_ui(kl, 1);
                fmpz_init(tot);
                for (l = 0; l < FLINT_MIN(20, fmpz_poly_length(b)); l++)
                {
                    fmpz_set(tot, il);
                    fmpz_add(tot, tot, jl);
                    fmpz_add(tot, tot, kl);
                    result = fmpz_equal(fmpz_poly_get_coeff_ptr(b, l), tot);
                    if (!result)
                    {
                        flint_printf("FAIL:\n\n");
                        flint_printf("%d %d %d %d\n", i, j, k, l);
                        fflush(stdout);
                        flint_abort();
                    }
                    fmpz_mul_si(il, il, i);
                    fmpz_mul_si(jl, jl, j);
                    fmpz_mul_si(kl, kl, k);
                }

                fmpz_poly_power_sums(b, a, 4);
                fmpz_poly_power_sums_to_poly(c, b);
                result = fmpz_poly_equal(a, c);
                if (!result)
                {
                    flint_printf("FAIL: newton series to poly\n\n");
                    fmpz_poly_print(a), flint_printf("\n\n");
                    fmpz_poly_print(b), flint_printf("\n\n");
                    fmpz_poly_print(c), flint_printf("\n\n");
                }

                fmpz_clear(il);
                fmpz_clear(jl);
                fmpz_clear(kl);
                fmpz_clear(tot);
                fmpz_poly_clear(a);
                fmpz_poly_clear(b);
                fmpz_poly_clear(c);
            }

    /* Check that the various implementations coincide and that */
    /* power_sums_to_poly gives back the original polynomial */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_init(e);
        fmpz_poly_init(f);
        fmpz_poly_randtest_not_zero(a, state, 1 + n_randint(state, 50), 100);
        fmpz_poly_set_coeff_ui(a, fmpz_poly_degree(a), 1);

        for (n = -1; n < 3 * fmpz_poly_degree(a);
             n += 1 + n_randint(state, fmpz_poly_degree(a)))
        {
            fmpz_poly_power_sums(b, a, n);
            fmpz_poly_power_sums_naive(c, a, n);

            result = fmpz_poly_equal(b, c);
            if (!result)
            {
                flint_printf
                    ("FAIL: equality power_sums, power_sums_naive\n");
                flint_printf("%ld", n);
                fmpz_poly_print(a), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            if (n >= 1)
            {
                /* use the formula PowerSums(p) = rev(poly') / rev(poly) */
                fmpz_poly_reverse(e, a, fmpz_poly_length(a));
                fmpz_poly_derivative(f, a);
                fmpz_poly_reverse(f, f, fmpz_poly_length(f));
                fmpz_poly_div_series(d, f, e, n);
                result = fmpz_poly_equal(b, d);
                if (!result)
                {
                    flint_printf("FAIL: equality with schoenhage formula\n");
                    flint_printf("%ld", n);
                    fmpz_poly_print(a), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_poly_power_sums(b, a, fmpz_poly_length(a));
        fmpz_poly_power_sums_to_poly(c, b);
        result = fmpz_poly_equal(a, c);
        if (!result)
        {
            flint_printf("FAIL: power_sums_to_poly\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
        fmpz_poly_clear(e);
        fmpz_poly_clear(f);
    }

    /* Check that the product of polynomials correspond to the sum of Power sums series */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_init(e);
        fmpz_poly_init(f);

        fmpz_poly_randtest_not_zero(a, state, 1 + n_randint(state, 5), 20);
        fmpz_poly_set_coeff_ui(a, fmpz_poly_degree(a), 1);
        fmpz_poly_randtest_not_zero(b, state, 1 + n_randint(state, 5), 20);
        fmpz_poly_set_coeff_ui(b, fmpz_poly_degree(b), 1);

        fmpz_poly_power_sums(c, a, 20);
        fmpz_poly_power_sums(d, b, 20);
        fmpz_poly_add(f, c, d);

        fmpz_poly_mul(e, a, b);
        fmpz_poly_power_sums(e, e, 20);

        result = fmpz_poly_equal(e, f);
        if (!result)
        {
            flint_printf("FAIL: PowerSums(p1 p2) = PowerSums(p1) + PowerSums(p2)\n");
            fmpz_poly_print(a), flint_printf("\n");
            fmpz_poly_print(b), flint_printf("\n");
            fmpz_poly_print(c), flint_printf("\n");
            fmpz_poly_print(d), flint_printf("\n");
            fmpz_poly_print(e), flint_printf("\n");
            fmpz_poly_print(f), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
        fmpz_poly_clear(e);
        fmpz_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
