/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

static void check_factorization(
    const fmpz_poly_factor_t fac,
    const fmpz_poly_t f,
    slong omega_lower_bound)
{
    slong i;
    slong omega;
    fmpz_t c;
    fmpz_poly_t h, t;

    fmpz_init(c);
    fmpz_poly_init(h);
    fmpz_poly_init(t);

    omega = 0;
    fmpz_poly_set_fmpz(h, &fac->c);
    for (i = 0; i < fac->num; i++)
    {
        omega += fac->exp[i];
        fmpz_poly_pow(t, fac->p + i, fac->exp[i]);
        fmpz_poly_mul(h, h, t);

        if (fac->p[i].length < 2)
        {
            flint_printf("FAIL\ncheck constant factor\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_content(c, fac->p + i);
        if (!fmpz_is_one(c) || fmpz_sgn(fac->p[i].coeffs + fac->p[i].length - 1) < 0)
        {
            flint_printf("FAIL\ncheck factor content\n");
            fflush(stdout);
            flint_abort();
        }
    }

    if (!fmpz_poly_equal(h, f))
    {
        flint_printf("FAIL\ncheck factorization matches\n");
        fflush(stdout);
        flint_abort();
    }

    if (omega < omega_lower_bound)
    {
        flint_printf("FAIL\ncheck omega\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_clear(c);
    fmpz_poly_clear(h);
    fmpz_poly_clear(t);
}

TEST_FUNCTION_START(fmpz_poly_factor_cubic, state)
{
    flint_bitcnt_t max_bits = 2000;
    slong i, tmul = 500;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        fmpz_poly_factor_t fac;

        fmpz_poly_init(g);
        fmpz_poly_init(f);
        fmpz_poly_factor_init(fac);

        /* (degree 1)^3 */
        do {
           fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
        } while (g->length != 2);
        fmpz_poly_pow(f, g, 3);
        fmpz_poly_factor(fac, f);
        check_factorization(fac, f, 3);

        /* (degree 1)*(degree 1)^2 */
        do {
           fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
        } while (g->length != 2);
        fmpz_poly_pow(f, g, 2);
        do {
           fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
        } while (g->length != 2);
        fmpz_poly_mul(f, f, g);
        fmpz_poly_factor(fac, f);
        check_factorization(fac, f, 3);

        /* (degree 1)*(degree 1)*(degree 1) */
        do {
           fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
        } while (g->length != 2);
        fmpz_poly_set(f, g);
        do {
           fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
        } while (g->length != 2);
        fmpz_poly_mul(f, f, g);
        do {
           fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
        } while (g->length != 2);
        fmpz_poly_mul(f, f, g);
        fmpz_poly_factor(fac, f);
        check_factorization(fac, f, 3);

        /* (degree 1)*(degree 2) */
        do {
           fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
        } while (g->length != 2);
        fmpz_poly_set(f, g);
        do {
           fmpz_poly_randtest(g, state, 3, n_randint(state, max_bits) + 2);
        } while (g->length != 3);
        fmpz_poly_mul(f, f, g);
        fmpz_poly_factor(fac, f);
        check_factorization(fac, f, 2);

        /* (degree 3) */
        do {
           fmpz_poly_randtest(g, state, 4, n_randint(state, max_bits) + 2);
        } while (g->length != 4);
        fmpz_poly_set(f, g);
        fmpz_poly_factor(fac, f);
        check_factorization(fac, f, 1);

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_factor_clear(fac);
    }

    TEST_FUNCTION_END(state);
}
