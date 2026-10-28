/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_resultant_discriminant, state)
{
    slong i, j;

    /* Check quadratic polynomial */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, d, d1;
        const char * vars[] = {"x","a","b","c"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(d, ctx);
        fmpz_mpoly_init(d1, ctx);

        fmpz_mpoly_set_str_pretty(f, "a^10*x^2 + b^100*x + c^100000000000000000000", vars, ctx);
        fmpz_mpoly_set_str_pretty(d1, "b^200 - 4*a^10*c^100000000000000000000", vars, ctx);
        if (!fmpz_mpoly_discriminant(d, f, 0, ctx))
        {
            flint_printf("FAIL: could not compute quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_mpoly_equal(d, d1, ctx))
        {
            flint_printf("FAIL: Check quadratic polynomial\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(d, ctx);
        fmpz_mpoly_clear(d1, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check univariate resultant */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, r;
        fmpz_poly_t au, bu;
        fmpz_t ru;
        ordering_t ord;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(r, ctx);
        fmpz_poly_init(au);
        fmpz_poly_init(bu);
        fmpz_init(ru);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 50) + 1;
        exp_bound2 = n_randint(state, 50) + 1;
        coeff_bits = n_randint(state, 4);
        fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
        fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);

        fmpz_mpoly_get_fmpz_poly(au, a, 0, ctx);
        fmpz_mpoly_get_fmpz_poly(bu, b, 0, ctx);

        fmpz_poly_resultant(ru, au, bu);

        if (!fmpz_mpoly_resultant(r, a, b, 0, ctx) ||
            !fmpz_mpoly_equal_fmpz(r, ru, ctx))
        {
            flint_printf("FAIL: Check univariate resultant \n");
            flint_printf("i: %wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_poly_clear(au);
        fmpz_poly_clear(bu);
        fmpz_clear(ru);
    }

    /* Check res(a*b,c) = res(a,c)*res(b,c) */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, ab, ra, rb, rab, p;
        ordering_t ord;
        slong nvars, len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 3) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(ab, ctx);
        fmpz_mpoly_init(ra, ctx);
        fmpz_mpoly_init(rb, ctx);
        fmpz_mpoly_init(rab, ctx);
        fmpz_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        len3 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 5) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        exp_bound3 = n_randint(state, 5) + 1;
        coeff_bits = n_randint(state, 3);
        fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
        fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);
        fmpz_mpoly_randtest_bound(c, state, len3, coeff_bits, exp_bound3, ctx);

        for (j = 0; j < nvars; j++)
        {
            fmpz_mpoly_mul(ab, a, b, ctx);

            if (!fmpz_mpoly_resultant(ra, a, c, j, ctx))
                continue;
            fmpz_mpoly_assert_canonical(ra, ctx);

            if (!fmpz_mpoly_resultant(rb, b, c, j, ctx))
                continue;
            fmpz_mpoly_assert_canonical(rb, ctx);

            if (!fmpz_mpoly_resultant(rab, ab, c, j, ctx))
                continue;
            fmpz_mpoly_assert_canonical(rab, ctx);

            fmpz_mpoly_mul(p, ra, rb, ctx);

            if (!fmpz_mpoly_equal(p,rab,ctx))
            {
                flint_printf("FAIL: Check res(a*b,c) = res(a,c)*res(b,c)\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(ab, ctx);
        fmpz_mpoly_clear(ra, ctx);
        fmpz_mpoly_clear(rb, ctx);
        fmpz_mpoly_clear(rab, ctx);
        fmpz_mpoly_clear(p, ctx);
    }

    /* Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, ab, r, da, db, dab, p;
        ordering_t ord;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 3) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(ab, ctx);
        fmpz_mpoly_init(da, ctx);
        fmpz_mpoly_init(db, ctx);
        fmpz_mpoly_init(dab, ctx);
        fmpz_mpoly_init(r, ctx);
        fmpz_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 5) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        coeff_bits = n_randint(state, 3);
        fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
        fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);

        for (j = 0; j < nvars; j++)
        {
            if (fmpz_mpoly_degree_si(a, j, ctx) < 1)
                continue;
            if (fmpz_mpoly_degree_si(b, j, ctx) < 1)
                continue;

            fmpz_mpoly_mul(ab, a, b, ctx);

            if (!fmpz_mpoly_resultant(r, a, b, j, ctx))
                continue;
            if (!fmpz_mpoly_discriminant(da, a, j, ctx))
                continue;
            if (!fmpz_mpoly_discriminant(db, b, j, ctx))
                continue;
            if (!fmpz_mpoly_discriminant(dab, ab, j, ctx))
                continue;

            fmpz_mpoly_mul(p, da, db, ctx);
            fmpz_mpoly_mul(p, p, r, ctx);
            fmpz_mpoly_mul(p, p, r, ctx);

            if (!fmpz_mpoly_equal(dab, p, ctx))
            {
                flint_printf("FAIL: Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(ab, ctx);
        fmpz_mpoly_clear(da, ctx);
        fmpz_mpoly_clear(db, ctx);
        fmpz_mpoly_clear(dab, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_mpoly_clear(p, ctx);
    }

    /* Check bivariate resultants, which take the dense algorithm at these
       sizes, against the univariate form, for inputs whose exponent fields
       are narrower than the resultant needs, inputs whose exponent fields
       take more than one word, and inputs with a common factor */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, r, s;
        fmpz_mpoly_univar_t au, bu;
        ulong bounds[2], exp[2];
        flint_bitcnt_t coeff_bits = 1 + n_randint(state, 40);
        slong var, dy, dx;
        int which = n_randint(state, 3);

        fmpz_mpoly_ctx_init(ctx, 2, mpoly_ordering_randtest(state));

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(r, ctx);
        fmpz_mpoly_init(s, ctx);
        fmpz_mpoly_univar_init(au, ctx);
        fmpz_mpoly_univar_init(bu, ctx);

        /* the resultant has degree about 2 dy dx > 255 in the other variable */
        /* the inputs have dy + 1 >= 7 coefficients in var, so that 
         *       (dy + 1)^2 >= 49 >= MULTIPOINT_CUTOFF
         * make sure to increase dy if MULTIPOINT_CUTOFF gets raised  */
        var = n_randint(state, 2);
        dy = 6 + n_randint(state, 3);
        dx = 22 + n_randint(state, 8);

        bounds[var] = dy + 1;
        bounds[1 - var] = dx + 1;
        exp[var] = dy;
        exp[1 - var] = dx;

        fmpz_mpoly_randtest_bounds(a, state, 15 + n_randint(state, 30), coeff_bits, bounds, ctx);
        fmpz_mpoly_randtest_bounds(b, state, 15 + n_randint(state, 30), coeff_bits, bounds, ctx);
        fmpz_mpoly_set_coeff_ui_ui(a, 1, exp, ctx);
        fmpz_mpoly_set_coeff_ui_ui(b, 1, exp, ctx);

        if (which == 0)
        {
            fmpz_mpoly_repack_bits_inplace(a, MPOLY_MIN_BITS, ctx);
            fmpz_mpoly_repack_bits_inplace(b, MPOLY_MIN_BITS, ctx);
        }
        else if (which == 1)
        {
            fmpz_mpoly_repack_bits_inplace(a, (2 + n_randint(state, 2)) * FLINT_BITS, ctx);
            fmpz_mpoly_repack_bits_inplace(b, (2 + n_randint(state, 2)) * FLINT_BITS, ctx);
        }
        else
        {
            bounds[var] = 2;
            bounds[1 - var] = 3;
            fmpz_mpoly_randtest_bounds(c, state, 4, coeff_bits, bounds, ctx);
            exp[var] = 1;
            exp[1 - var] = 0;
            fmpz_mpoly_set_coeff_ui_ui(c, 1, exp, ctx);
            fmpz_mpoly_mul(a, a, c, ctx);
            fmpz_mpoly_mul(b, b, c, ctx);
        }

        fmpz_mpoly_to_univar(au, a, var, ctx);
        fmpz_mpoly_to_univar(bu, b, var, ctx);

        if (fmpz_mpoly_univar_resultant(s, au, bu, ctx))
        {
            if (!fmpz_mpoly_resultant(r, a, b, var, ctx))
            {
                flint_printf("FAIL: could not compute bivariate resultant\n");
                flint_printf("i: %wd  which: %d\n", i, which);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mpoly_assert_canonical(r, ctx);

            if (!fmpz_mpoly_equal(r, s, ctx) ||
                (which == 2 && !fmpz_mpoly_is_zero(r, ctx)))
            {
                flint_printf("FAIL: Check bivariate resultant\n");
                flint_printf("i: %wd  which: %d\n", i, which);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_univar_clear(au, ctx);
        fmpz_mpoly_univar_clear(bu, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_mpoly_clear(s, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
