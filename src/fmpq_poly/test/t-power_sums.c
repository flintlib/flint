/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_power_sums, state)
{
    int i, j, k, l, den, result;
    fmpz_t il, jl, kl, dl;
    fmpq_t tot;
    fmpq_t tmp;
    fmpq_poly_t a, b, c, d;
    fmpz_poly_t az, bz;

    /* Check that it is valid in degree 3 with rational roots, ie */
    /* for polynomials of the form (dx-i)(dx-j)(dx-k)             */
    for (den = 1; den < 6; den++)
        for (i = -4; i < 4; i++)
            for (j = i; j < 4; j++)
                for (k = j; k < 4; k++)
                {
                    fmpq_poly_init(a);
                    fmpq_poly_init(b);
                    fmpq_poly_init(c);

                    fmpq_poly_set_coeff_si(a, 0, -i * j * k);
                    fmpq_poly_set_coeff_si(a, 1,
                                           den * (i * j + i * k + j * k));
                    fmpq_poly_set_coeff_si(a, 2, den * den * (-i - j - k));
                    fmpq_poly_set_coeff_si(a, 3, den * den * den);

                    fmpq_poly_power_sums(b, a, 20);

                    fmpz_init(il);
                    fmpz_init(jl);
                    fmpz_init(kl);
                    fmpz_init(dl);
                    fmpq_init(tot);
                    fmpq_init(tmp);
                    fmpz_one(il);
                    fmpz_one(jl);
                    fmpz_one(kl);
                    fmpz_one(dl);
                    for (l = 0; l < FLINT_MIN(20, fmpq_poly_length(b)); l++)
                    {
                        fmpq_zero(tot);
                        fmpq_add_fmpz(tot, tot, il);
                        fmpq_add_fmpz(tot, tot, jl);
                        fmpq_add_fmpz(tot, tot, kl);
                        fmpq_div_fmpz(tot, tot, dl);

                        fmpq_poly_get_coeff_fmpq(tmp, b, l);
                        result = fmpq_equal(tmp, tot);
                        if (!result)
                        {
                            flint_printf("FAIL: power sums rational root\n");
                            flint_printf("%d %d %d %d %d\n", i, j, k, den, l);
                            fflush(stdout);
                            flint_abort();
                        }
                        fmpz_mul_si(il, il, i);
                        fmpz_mul_si(jl, jl, j);
                        fmpz_mul_si(kl, kl, k);
                        fmpz_mul_si(dl, dl, den);
                    }

                    fmpq_poly_power_sums(b, a, 4);
                    fmpq_poly_power_sums_to_poly(c, b);
                    fmpq_poly_make_monic(a, a);
                    result = fmpq_poly_equal(a, c);
                    if (!result)
                    {
                        flint_printf("FAIL: power sums to poly\n");
                        fmpq_poly_print(a), flint_printf("\n\n");
                        fmpq_poly_print(b), flint_printf("\n\n");
                        fflush(stdout);
                        flint_abort();
                    }

                    fmpz_clear(il);
                    fmpz_clear(jl);
                    fmpz_clear(kl);
                    fmpq_clear(tot);
                    fmpq_clear(tmp);
                    fmpq_poly_clear(a);
                    fmpq_poly_clear(b);
                    fmpq_poly_clear(c);
                }

    /* Check that going back and forth between the power sums representation gives the identity */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_init(az);
        fmpz_poly_init(bz);
        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);

        fmpz_poly_randtest_not_zero(az, state, 1 + n_randint(state, 15), 30);
        fmpq_poly_set_fmpz_poly(a, az);

        fmpq_poly_power_sums(b, a, 20);
        fmpq_poly_power_sums_to_fmpz_poly(bz, b);
        fmpq_poly_power_sums_to_poly(c, b);

        fmpq_poly_make_monic(a, a);
        _fmpz_poly_primitive_part(az->coeffs, az->coeffs, az->length);

        result = fmpq_poly_equal(a, c) && fmpz_poly_equal(bz, az);
        if (!result)
        {
            flint_printf("FAIL: power sums - power sums to poly\n");
            fmpz_poly_print(az), flint_printf("\n\n");
            fmpz_poly_print(bz), flint_printf("\n\n");
            fmpq_poly_print(a), flint_printf("\n\n");
            fmpq_poly_print(b), flint_printf("\n\n");
            fmpq_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpz_poly_clear(az);
        fmpz_poly_clear(bz);
    }

    /* Check that the product of polynomials correspond to the sum of Power sums series */
    /* (and aliasing of fmpq_poly_power_sums)                                           */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(d);

        fmpq_poly_randtest_not_zero(a, state, 1 + n_randint(state, 10), 30);
        fmpq_poly_randtest_not_zero(b, state, 1 + n_randint(state, 10), 30);

        fmpq_poly_mul(c, a, b);
        fmpq_poly_power_sums(c, c, 20);

        fmpq_poly_power_sums(a, a, 20);
        fmpq_poly_power_sums(b, b, 20);
        fmpq_poly_add(d, a, b);

        result = fmpq_poly_equal(c, d);
        if (!result)
        {
            flint_printf
                ("FAIL: PowerSums(p1 p2) = PowerSums(p1) + PowerSums(p2)\n");
            fmpq_poly_print(a), flint_printf("\n");
            fmpq_poly_print(b), flint_printf("\n");
            fmpq_poly_print(c), flint_printf("\n");
            fmpq_poly_print(d), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
