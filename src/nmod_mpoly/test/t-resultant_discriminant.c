/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "nmod_poly.h"
#include "mpoly.h"
#include "nmod_mpoly.h"

/* keeps the terms of A of weighted degree w[0] e[0] + ... at most d */
static void
_nmod_mpoly_truncate_weighted(nmod_mpoly_t A, const ulong * w, ulong d,
                              const nmod_mpoly_ctx_t ctx)
{
    slong i, k, len = 0, nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    ulong * e = (ulong *) flint_malloc(nvars*sizeof(ulong));
    ulong s;

    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ui(e, A->exps + N*i, A->bits, ctx->minfo);

        s = 0;
        for (k = 0; k < nvars; k++)
            s += w[k]*e[k];

        if (s > d)
            continue;

        A->coeffs[len] = A->coeffs[i];
        mpoly_monomial_set(A->exps + N*len, A->exps + N*i, N);
        len++;
    }

    _nmod_mpoly_set_length(A, len, ctx);

    flint_free(e);
}

TEST_FUNCTION_START(nmod_mpoly_resultant_discriminant, state)
{
    slong i, j;

    /* Check quadratic polynomial */
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, d, d1;
        const char * vars[] = {"x","a","b","c"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_DEGLEX, 100003);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(d, ctx);
        nmod_mpoly_init(d1, ctx);

        nmod_mpoly_set_str_pretty(f, "a^10*x^2 + b^100*x + c^100000000000000000000", vars, ctx);
        nmod_mpoly_set_str_pretty(d1, "b^200 - 4*a^10*c^100000000000000000000", vars, ctx);
        if (!nmod_mpoly_discriminant(d, f, 0, ctx))
        {
            flint_printf("FAIL: could not compute quadratic discriminant\n");
            fflush(stdout);
            flint_abort();
        }

        if (!nmod_mpoly_equal(d, d1, ctx))
        {
            flint_printf("FAIL: Check quadratic polynomial\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(d, ctx);
        nmod_mpoly_clear(d1, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check univariate resultant */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, r;
        nmod_poly_t au, bu;
        ulong ru;
        slong len1, len2, exp_bound1, exp_bound2;
        ulong modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, 1, mpoly_ordering_randtest(state), modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(r, ctx);
        nmod_poly_init(au, nmod_mpoly_ctx_modulus(ctx));
        nmod_poly_init(bu, nmod_mpoly_ctx_modulus(ctx));

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 50) + 1;
        exp_bound2 = n_randint(state, 50) + 1;
        nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        nmod_mpoly_get_nmod_poly(au, a, 0, ctx);
        nmod_mpoly_get_nmod_poly(bu, b, 0, ctx);

        ru = nmod_poly_resultant(au, bu);

        if (!nmod_mpoly_resultant(r, a, b, 0, ctx) ||
            !nmod_mpoly_equal_ui(r, ru, ctx))
        {
            flint_printf("FAIL: Check univariate resultant \n");
            flint_printf("i: %wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_poly_clear(au);
        nmod_poly_clear(bu);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check res(a*b,c) = res(a,c)*res(b,c) */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, c, ab, ra, rb, rab, p;
        slong len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        ulong modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 3, modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(c, ctx);
        nmod_mpoly_init(ab, ctx);
        nmod_mpoly_init(ra, ctx);
        nmod_mpoly_init(rb, ctx);
        nmod_mpoly_init(rab, ctx);
        nmod_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        len3 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 5) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        exp_bound3 = n_randint(state, 5) + 1;
        nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);
        nmod_mpoly_randtest_bound(c, state, len3, exp_bound3, ctx);

        for (j = 0; j < nmod_mpoly_ctx_nvars(ctx); j++)
        {
            nmod_mpoly_mul(ab, a, b, ctx);

            if (!nmod_mpoly_resultant(ra, a, c, j, ctx))
                continue;
            nmod_mpoly_assert_canonical(ra, ctx);

            if (!nmod_mpoly_resultant(rb, b, c, j, ctx))
                continue;
            nmod_mpoly_assert_canonical(rb, ctx);

            if (!nmod_mpoly_resultant(rab, ab, c, j, ctx))
                continue;
            nmod_mpoly_assert_canonical(rab, ctx);

            nmod_mpoly_mul(p, ra, rb, ctx);

            if (!nmod_mpoly_equal(p,rab,ctx))
            {
                flint_printf("FAIL: Check res(a*b,c) = res(a,c)*res(b,c)\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(c, ctx);
        nmod_mpoly_clear(ab, ctx);
        nmod_mpoly_clear(ra, ctx);
        nmod_mpoly_clear(rb, ctx);
        nmod_mpoly_clear(rab, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, ab, r, da, db, dab, p;
        slong len1, len2, exp_bound1, exp_bound2;
        ulong modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 3, modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ab, ctx);
        nmod_mpoly_init(da, ctx);
        nmod_mpoly_init(db, ctx);
        nmod_mpoly_init(dab, ctx);
        nmod_mpoly_init(r, ctx);
        nmod_mpoly_init(p, ctx);

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        exp_bound1 = n_randint(state, 4) + 1;
        exp_bound2 = n_randint(state, 5) + 1;
        nmod_mpoly_randtest_bound(a, state, len1, exp_bound1, ctx);
        nmod_mpoly_randtest_bound(b, state, len2, exp_bound2, ctx);

        for (j = 0; j < nmod_mpoly_ctx_nvars(ctx); j++)
        {
            if (nmod_mpoly_degree_si(a, j, ctx) < 1)
                continue;
            if (nmod_mpoly_degree_si(b, j, ctx) < 1)
                continue;

            nmod_mpoly_mul(ab, a, b, ctx);

            if (!nmod_mpoly_resultant(r, a, b, j, ctx))
                continue;
            if (!nmod_mpoly_discriminant(da, a, j, ctx))
                continue;
            if (!nmod_mpoly_discriminant(db, b, j, ctx))
                continue;
            if (!nmod_mpoly_discriminant(dab, ab, j, ctx))
                continue;

            nmod_mpoly_mul(p, da, db, ctx);
            nmod_mpoly_mul(p, p, r, ctx);
            nmod_mpoly_mul(p, p, r, ctx);

            if (!nmod_mpoly_equal(dab, p, ctx))
            {
                flint_printf("FAIL: Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2\n");

                flint_printf("a: "); nmod_mpoly_print_pretty(a, NULL, ctx); flint_printf("\n");
                flint_printf("b: "); nmod_mpoly_print_pretty(b, NULL, ctx); flint_printf("\n");

                flint_printf("disc(a*b): "); nmod_mpoly_print_pretty(dab, NULL, ctx); flint_printf("\n");
                flint_printf("disc(a): "); nmod_mpoly_print_pretty(da, NULL, ctx); flint_printf("\n");
                flint_printf("disc(b): "); nmod_mpoly_print_pretty(db, NULL, ctx); flint_printf("\n");
                flint_printf("res(a, b): "); nmod_mpoly_print_pretty(r, NULL, ctx); flint_printf("\n");

                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(ab, ctx);
        nmod_mpoly_clear(da, ctx);
        nmod_mpoly_clear(db, ctx);
        nmod_mpoly_clear(dab, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check bivariate resultants, which take the dense algorithm at these
       sizes, against the univariate form, for inputs whose exponent fields
       are narrower than the resultant needs, inputs whose exponent fields
       take more than one word, and inputs with a common factor */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, c, r, s;
        nmod_mpoly_univar_t au, bu;
        ulong modulus, bounds[2], exp[2];
        slong var, dy, dx;
        int which = n_randint(state, 3);

        /* large enough for the number of points of the dense algorithm */
        modulus = n_randint(state, FLINT_BITS - 20) + 20;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init(ctx, 2, mpoly_ordering_randtest(state), modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(c, ctx);
        nmod_mpoly_init(r, ctx);
        nmod_mpoly_init(s, ctx);
        nmod_mpoly_univar_init(au, ctx);
        nmod_mpoly_univar_init(bu, ctx);

        /* the resultant has degree about 2 dy dx > 255 in the other variable */
        var = n_randint(state, 2);
        dy = 6 + n_randint(state, 4);
        dx = 22 + n_randint(state, 10);

        bounds[var] = dy + 1;
        bounds[1 - var] = dx + 1;
        exp[var] = dy;
        exp[1 - var] = dx;

        nmod_mpoly_randtest_bounds(a, state, 20 + n_randint(state, 60), bounds, ctx);
        nmod_mpoly_randtest_bounds(b, state, 20 + n_randint(state, 60), bounds, ctx);
        nmod_mpoly_set_coeff_ui_ui(a, 1, exp, ctx);
        nmod_mpoly_set_coeff_ui_ui(b, 1, exp, ctx);

        if (which == 0)
        {
            nmod_mpoly_repack_bits_inplace(a, MPOLY_MIN_BITS, ctx);
            nmod_mpoly_repack_bits_inplace(b, MPOLY_MIN_BITS, ctx);
        }
        else if (which == 1)
        {
            nmod_mpoly_repack_bits_inplace(a, (2 + n_randint(state, 2)) * FLINT_BITS, ctx);
            nmod_mpoly_repack_bits_inplace(b, (2 + n_randint(state, 2)) * FLINT_BITS, ctx);
        }
        else
        {
            bounds[var] = 2;
            bounds[1 - var] = 3;
            nmod_mpoly_randtest_bounds(c, state, 4, bounds, ctx);
            exp[var] = 1;
            exp[1 - var] = 0;
            nmod_mpoly_set_coeff_ui_ui(c, 1, exp, ctx);
            nmod_mpoly_mul(a, a, c, ctx);
            nmod_mpoly_mul(b, b, c, ctx);
        }

        nmod_mpoly_to_univar(au, a, var, ctx);
        nmod_mpoly_to_univar(bu, b, var, ctx);

        if (nmod_mpoly_univar_resultant(s, au, bu, ctx))
        {
            if (!nmod_mpoly_resultant(r, a, b, var, ctx))
            {
                flint_printf("FAIL: could not compute bivariate resultant\n");
                flint_printf("i: %wd  which: %d\n", i, which);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_assert_canonical(r, ctx);

            if (!nmod_mpoly_equal(r, s, ctx) ||
                (which == 2 && !nmod_mpoly_is_zero(r, ctx)))
            {
                flint_printf("FAIL: Check bivariate resultant\n");
                flint_printf("i: %wd  which: %d\n", i, which);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_univar_clear(au, ctx);
        nmod_mpoly_univar_clear(bu, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(c, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_mpoly_clear(s, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check resultants in several variables, which go through Kronecker
       substitution into the dense bivariate algorithm, against the
       univariate form: with exponents that deflate, monomial factors,
       variables absent from both inputs, common factors, narrow exponent
       fields, aliasing, and fields too small for the dense algorithm */
    {
        slong save_threads = flint_get_num_threads();

        for (i = 0; i < 40 * flint_test_multiplier(); i++)
        {
            nmod_mpoly_ctx_t ctx;
            nmod_mpoly_t a, b, c, r, s;
            nmod_mpoly_univar_t au, bu;
            ulong modulus, * bounds, * exp;
            fmpz * shift, * stride;
            slong nvars, var, j, k;
            int common = (n_randint(state, 4) == 0);

            /* now and then too small for the number of evaluation points */
            if (n_randint(state, 4) == 0)
                modulus = n_randint(state, FLINT_BITS - 1) + 1;
            else
                modulus = n_randint(state, FLINT_BITS - 20) + 20;
            modulus = n_randbits(state, modulus);
            modulus = n_nextprime(modulus, 1);

            nvars = 2 + n_randint(state, 3);
            nmod_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state), modulus);

            nmod_mpoly_init(a, ctx);
            nmod_mpoly_init(b, ctx);
            nmod_mpoly_init(c, ctx);
            nmod_mpoly_init(r, ctx);
            nmod_mpoly_init(s, ctx);
            nmod_mpoly_univar_init(au, ctx);
            nmod_mpoly_univar_init(bu, ctx);

            bounds = (ulong *) flint_malloc(2*nvars*sizeof(ulong));
            exp = bounds + nvars;
            shift = _fmpz_vec_init(nvars);
            stride = _fmpz_vec_init(nvars);

            flint_set_num_threads(1 + n_randint(state, 4));

            var = n_randint(state, nvars);
            for (k = 0; k < nvars; k++)
            {
                if (k == var)
                {
                    bounds[k] = 2 + n_randint(state, 4);
                    fmpz_one(stride + k);
                }
                else
                {
                    bounds[k] = (n_randint(state, 4) == 0) ? 1 : 2 + n_randint(state, 3);
                    /* a common factor undoes the deflation, which would
                       only leave larger inputs */
                    fmpz_set_ui(stride + k, common ? 1 : 1 + n_randint(state, 3));
                    fmpz_set_ui(shift + k, n_randint(state, 3));
                }
            }

            nmod_mpoly_randtest_bounds(a, state, 1 + n_randint(state, 40), bounds, ctx);
            nmod_mpoly_randtest_bounds(b, state, 1 + n_randint(state, 40), bounds, ctx);

            /* half of the time only the terms of bounded weighted degree, for
               which bounds on the degrees of the resultant of weight above 0
               in y are the smaller ones */
            if (n_randint(state, 2))
            {
                ulong d;

                for (j = 0; j < 2; j++)
                {
                    d = 0;
                    for (k = 0; k < nvars; k++)
                    {
                        exp[k] = 1 + n_randint(state, 3);
                        d += exp[k]*(bounds[k] - 1);
                    }
                    d = d/2 + n_randint(state, d/2 + 1);
                    _nmod_mpoly_truncate_weighted(j == 0 ? a : b, exp, d, ctx);
                }
            }

            /* exponents e -> shift + stride e in the variables other than var,
               with a different monomial factor for b */
            nmod_mpoly_inflate(a, a, shift, stride, ctx);
            k = (var + 1) % nvars;
            fmpz_add_ui(shift + k, shift + k, 1);
            nmod_mpoly_inflate(b, b, shift, stride, ctx);

            if (common)
            {
                for (k = 0; k < nvars; k++)
                {
                    bounds[k] = 2;
                    exp[k] = (k == var);
                }
                nmod_mpoly_randtest_bounds(c, state, 3, bounds, ctx);
                nmod_mpoly_set_coeff_ui_ui(c, 1, exp, ctx);
                nmod_mpoly_mul(a, a, c, ctx);
                nmod_mpoly_mul(b, b, c, ctx);
            }

            if (n_randint(state, 2))
            {
                nmod_mpoly_repack_bits_inplace(a, MPOLY_MIN_BITS, ctx);
                nmod_mpoly_repack_bits_inplace(b, MPOLY_MIN_BITS, ctx);
            }

            nmod_mpoly_to_univar(au, a, var, ctx);
            nmod_mpoly_to_univar(bu, b, var, ctx);

            if (nmod_mpoly_univar_resultant(s, au, bu, ctx))
            {
                if (nmod_mpoly_resultant_kronecker(r, a, b, var, ctx))
                {
                    nmod_mpoly_assert_canonical(r, ctx);

                    if (!nmod_mpoly_equal(r, s, ctx))
                    {
                        flint_printf("FAIL: Check Kronecker resultant\n");
                        flint_printf("i: %wd  nvars: %wd  var: %wd\n", i, nvars, var);
                        fflush(stdout);
                        flint_abort();
                    }
                }

                nmod_mpoly_set(r, a, ctx);
                if (!nmod_mpoly_resultant(r, r, b, var, ctx))
                {
                    flint_printf("FAIL: could not compute multivariate resultant\n");
                    flint_printf("i: %wd  nvars: %wd  var: %wd\n", i, nvars, var);
                    fflush(stdout);
                    flint_abort();
                }

                nmod_mpoly_assert_canonical(r, ctx);

                if (!nmod_mpoly_equal(r, s, ctx) ||
                    (common && !nmod_mpoly_is_zero(r, ctx)))
                {
                    flint_printf("FAIL: Check multivariate resultant\n");
                    flint_printf("i: %wd  nvars: %wd  var: %wd\n", i, nvars, var);
                    fflush(stdout);
                    flint_abort();
                }
            }

            flint_free(bounds);
            _fmpz_vec_clear(shift, nvars);
            _fmpz_vec_clear(stride, nvars);
            nmod_mpoly_univar_clear(au, ctx);
            nmod_mpoly_univar_clear(bu, ctx);
            nmod_mpoly_clear(a, ctx);
            nmod_mpoly_clear(b, ctx);
            nmod_mpoly_clear(c, ctx);
            nmod_mpoly_clear(r, ctx);
            nmod_mpoly_clear(s, ctx);
            nmod_mpoly_ctx_clear(ctx);
        }

        flint_set_num_threads(save_threads);
    }

    TEST_FUNCTION_END(state);
}
