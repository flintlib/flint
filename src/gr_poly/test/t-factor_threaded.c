/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"

/* The parallel code paths are only used for inputs longer than
   gr_poly_factor_threaded_cutoff; the cutoff is lowered here so that
   they can be tested with small inputs. */
FLINT_DLL extern slong gr_poly_factor_threaded_cutoff;

#define MIN_LEN 24

/* Product of nfac distinct random monic irreducible polynomials of
   degree d (the factorization algorithms require squarefree input). */
static int
_random_equal_deg(gr_poly_t res, slong d, slong nfac, flint_rand_t state, gr_ctx_t ctx)
{
    gr_poly_t A, G;
    slong i;
    int status = GR_SUCCESS;

    gr_poly_init(A, ctx);
    gr_poly_init(G, ctx);
    status |= gr_poly_one(res, ctx);

    for (i = 0; i < nfac; i++)
    {
        do {
            status |= gr_poly_randtest(A, state, d, ctx);
            status |= gr_poly_set_coeff_ui(A, d, 1, ctx);

            if (gr_poly_is_irreducible(A, ctx) != T_TRUE)
                continue;

            status |= gr_poly_gcd(G, res, A, ctx);
        } while (status != GR_SUCCESS || gr_poly_is_irreducible(A, ctx) != T_TRUE
                    || G->length != 1);

        status |= gr_poly_mul(res, res, A, ctx);
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(G, ctx);
    return status;
}

TEST_FUNCTION_START(gr_poly_factor_threaded, state)
{
    slong iter;
    slong save_cutoff = gr_poly_factor_threaded_cutoff;

    gr_poly_factor_threaded_cutoff = 8;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t P, Q;
        gr_poly_vec_t fac1, fac2, dd1, dd2;
        fmpz_vec_t exp1, exp2, degs1, degs2;
        gr_ptr c1, c2;
        slong i, j, d, nfac, nthreads;
        int status = GR_SUCCESS;
        int found;

        /* small fields, so that the degrees needed to reach the
           threaded code paths stay cheap */
        /* The inputs below are products of distinct irreducible
           polynomials; when many of the same (small) degree are needed,
           the field has to be large enough to contain them. */
        if (iter % 2 == 0)
            gr_ctx_init_nmod(ctx, n_randprime(state, 13, 1));
        else
            switch (n_randint(state, 3))
            {
                case 0: gr_ctx_init_nmod(ctx, 2); break;
                case 1: gr_ctx_init_nmod(ctx, n_randprime(state, 5, 1)); break;
                default: gr_ctx_init_fq_nmod(ctx, n_randprime(state, 4, 1), 2, "a"); break;
            }

        GR_IGNORE(gr_ctx_set_is_field(ctx, T_TRUE));   /* fq contexts already know */

        gr_poly_init(P, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_vec_init(fac1, 0, ctx);
        gr_poly_vec_init(fac2, 0, ctx);
        gr_poly_vec_init(dd1, 0, ctx);
        gr_poly_vec_init(dd2, 0, ctx);
        fmpz_vec_init(exp1, 0);
        fmpz_vec_init(exp2, 0);
        fmpz_vec_init(degs1, 0);
        fmpz_vec_init(degs2, 0);
        c1 = gr_heap_init(ctx);
        c2 = gr_heap_init(ctx);

        if (iter % 4 == 0)
        {
            /* many small factors, so that the pieces resulting from the
               first splits are still long enough to be processed in
               parallel; d >= 16 also selects the doubling algorithm for
               the trace */
            /* d cycles through 1 (no trace) and values above 16, where
               the trace uses the doubling algorithm, with different bit
               patterns (even, one bit, several bits) */
            {
                static const slong ds[4] = { 1, 16, 17, 19 };
                d = ds[(iter / 4) % 4];
            }
            nfac = (d == 1) ? 2 * MIN_LEN + n_randint(state, 4) : 2 + n_randint(state, 2);
            status |= _random_equal_deg(P, d, nfac, state, ctx);
        }
        else if (iter % 2 == 0)
        {
            /* several factors of the same degree: exercises the
               parallel equal degree factorization */
            d = 1 + n_randint(state, 3);
            nfac = MIN_LEN / d + 1 + n_randint(state, 4);
            status |= _random_equal_deg(P, d, nfac, state, ctx);
        }
        else
        {
            /* factors of many different degrees: exercises the
               parallel distinct degree factorization. Half of the time
               only two factors of large degree, so that the distinct
               degree loop runs over several blocks of giant steps. */
            gr_poly_t A;
            gr_poly_init(A, ctx);
            status |= gr_poly_one(P, ctx);

            if (iter % 8 == 1)
            {
                /* two factors of large degree: the distinct degree loop
                   then runs over several blocks of giant steps */
                for (d = MIN_LEN; d <= MIN_LEN + 1; d++)
                {
                    status |= _random_equal_deg(A, d, 1, state, ctx);
                    status |= gr_poly_mul(P, P, A, ctx);
                }
            }
            else
            {
                for (d = 1; P->length - 1 < MIN_LEN; d++)
                {
                    status |= _random_equal_deg(A, d, 1, state, ctx);
                    status |= gr_poly_mul(P, P, A, ctx);
                }
            }

            gr_poly_clear(A, ctx);
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        nthreads = 2 + n_randint(state, 3);

        /* distinct degree factorization, serial and threaded */
        flint_set_num_threads(1);
        status |= gr_poly_factor_distinct_deg(dd1, degs1, P, ctx);
        flint_set_num_threads(nthreads);
        status |= gr_poly_factor_distinct_deg(dd2, degs2, P, ctx);

        if (status == GR_SUCCESS && dd1->length != dd2->length)
        {
            flint_printf("FAIL (distinct deg: number of pieces)\n\n");
            gr_ctx_println(ctx);
            flint_printf("threads = %wd, deg = %wd: %wd vs %wd\n",
                nthreads, P->length - 1, dd1->length, dd2->length);
            flint_abort();
        }

        for (i = 0; i < dd1->length && status == GR_SUCCESS; i++)
        {
            found = 0;
            for (j = 0; j < dd2->length; j++)
                if (gr_poly_equal(dd1->entries + i, dd2->entries + j, ctx) == T_TRUE
                        && fmpz_equal(degs1->entries + i, degs2->entries + j))
                    found = 1;

            if (!found)
            {
                flint_printf("FAIL (distinct deg: piece not found)\n\n");
                gr_ctx_println(ctx);
                flint_printf("threads = %wd, piece %wd of degree %wd\n",
                    nthreads, i, dd1->entries[i].length - 1);
                flint_abort();
            }
        }

        /* complete factorization, serial and threaded */
        flint_set_num_threads(1);
        status |= gr_poly_factor_finite_field(c1, fac1, exp1, P, ctx);
        flint_set_num_threads(nthreads);
        status |= gr_poly_factor_finite_field(c2, fac2, exp2, P, ctx);

        if (status == GR_SUCCESS)
        {
            if (fac1->length != fac2->length || gr_equal(c1, c2, ctx) == T_FALSE)
            {
                flint_printf("FAIL (factor: number of factors)\n\n");
                gr_ctx_println(ctx);
                flint_printf("threads = %wd, deg = %wd: %wd vs %wd\n",
                    nthreads, P->length - 1, fac1->length, fac2->length);
                flint_abort();
            }

            /* the threaded factorization multiplies out to the input */
            status |= gr_poly_set_scalar(Q, c2, ctx);
            for (i = 0; i < fac2->length; i++)
            {
                gr_poly_t A;
                gr_poly_init(A, ctx);
                status |= gr_poly_pow_ui(A, fac2->entries + i, fmpz_get_ui(exp2->entries + i), ctx);
                status |= gr_poly_mul(Q, Q, A, ctx);
                gr_poly_clear(A, ctx);

                if (gr_poly_is_irreducible(fac2->entries + i, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (factor: reducible factor)\n\n");
                    gr_ctx_println(ctx);
                    flint_abort();
                }
            }

            if (status == GR_SUCCESS && gr_poly_equal(P, Q, ctx) == T_FALSE)
            {
                flint_printf("FAIL (factor: product)\n\n");
                gr_ctx_println(ctx);
                flint_printf("threads = %wd, deg = %wd\n", nthreads, P->length - 1);
                flint_abort();
            }
        }

cleanup:
        flint_set_num_threads(1);

        gr_poly_clear(P, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_vec_clear(fac1, ctx);
        gr_poly_vec_clear(fac2, ctx);
        gr_poly_vec_clear(dd1, ctx);
        gr_poly_vec_clear(dd2, ctx);
        fmpz_vec_clear(exp1);
        fmpz_vec_clear(exp2);
        fmpz_vec_clear(degs1);
        fmpz_vec_clear(degs2);
        gr_heap_clear(c1, ctx);
        gr_heap_clear(c2, ctx);
        gr_ctx_clear(ctx);
    }

    gr_poly_factor_threaded_cutoff = save_cutoff;

    TEST_FUNCTION_END(state);
}
