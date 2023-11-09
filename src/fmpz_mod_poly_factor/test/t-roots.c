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
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"

/* Defined in t-roots.c and t-roots_factored.c */
#define test_poly test_poly_roots
void test_poly(
    fmpz_mod_poly_factor_t roots,
    const fmpz_mod_poly_t f,
    int want_mult,
    const fmpz_mod_ctx_t ctx)
{
    slong i, multiplicity;
    fmpz_mod_poly_t q, qt, r;

    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_init(qt, ctx);
    fmpz_mod_poly_init(r, ctx);
    fmpz_mod_poly_set(q, f, ctx);

    fmpz_mod_poly_roots(roots, f, want_mult, ctx);

    for (i = 0; i < roots->num; i++)
    {
        if (fmpz_mod_poly_degree(roots->poly + i, ctx) != 1)
        {
            flint_printf("FAILED:\ncheck root is linear\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_is_one(roots->poly[i].coeffs + 1))
        {
            flint_printf("FAILED:\ncheck root is monic\n");
            fflush(stdout);
            flint_abort();
        }

        multiplicity = 0;
        while (fmpz_mod_poly_divrem(qt, r, q, roots->poly + i, ctx),
               fmpz_mod_poly_is_zero(r, ctx))
        {
            fmpz_mod_poly_swap(q, qt, ctx);
            multiplicity++;
        }

        if (multiplicity <= 0)
        {
            flint_printf("FAILED:\ncheck root is a root\n");
            fflush(stdout);
            flint_abort();
        }

        if (roots->exp[i] != (want_mult ? multiplicity : 1))
        {
            flint_printf("FAILED:\ncheck root multiplicity\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_mod_poly_roots(roots, q, want_mult, ctx);
    if (roots->num > 0)
    {
        flint_printf("FAILED:\ncheck missing roots\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_poly_clear(qt, ctx);
    fmpz_mod_poly_clear(r, ctx);
}

TEST_FUNCTION_START(fmpz_mod_poly_factor_roots, state)
{
    slong i, j, k, l;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t f;
        fmpz_mod_poly_factor_t r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 100);

        fmpz_nextprime(p, p, 1);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_factor_init(r, ctx);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mod_poly_randtest(f, state, n_randint(state, 20) + 1, ctx);
            } while (fmpz_mod_poly_is_zero(f, ctx));

            for (k = 0; k < 5; k++)
            {
                fmpz_mod_poly_t ff;
                fmpz_mod_poly_init(ff, ctx);
                fmpz_mod_poly_fit_length(ff, 2, ctx);
                fmpz_one(ff->coeffs + 1);
                fmpz_randm(ff->coeffs + 0, state, p);
                ff->length = 2;
                for (l = 1 + n_randint(state, 5); l > 0; l--)
                    fmpz_mod_poly_mul(f, f, ff, ctx);
                fmpz_mod_poly_clear(ff, ctx);
            }

            if (n_randint(state, 2))
            {
                test_poly(r, f, 1, ctx);
                test_poly(r, f, 0, ctx);
            }
            else
            {
                test_poly(r, f, 0, ctx);
                test_poly(r, f, 1, ctx);
            }
        }

        fmpz_mod_poly_factor_clear(r, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
#undef test_poly
