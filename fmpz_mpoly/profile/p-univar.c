/*
    Copyright 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "fmpz_mpoly.h"


void run_cmp(
    fmpz_mpoly_struct * mpolys,
    fmpz_poly_struct * polys,
    fmpz_mpoly_struct * mres,
    fmpz_poly_struct * res,
    slong npolys,
    const fmpz_mpoly_ctx_t ctx,
    void (*f)(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t),
    void (*g)(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t))
{
    slong i, j, k;
    slong nreps = 100;
    timeit_t timer;

    for (j = 0; j < 2; j++)
    {
        timeit_start(timer);
        for (k = 0; k < nreps; k++)
        for (i = 1; i < npolys; i++)
            f(mres + i, mpolys + i, mpolys + i - 1, ctx);
        timeit_stop(timer);
        flint_printf("mpoly: %wd  ", timer->wall);
        fflush(stdout);

        timeit_start(timer);
        for (k = 0; k < nreps; k++)
        for (i = 1; i < npolys; i++)
            g(res + i, polys + i, polys + i - 1);
        timeit_stop(timer);
        flint_printf("poly: %wd\n", timer->wall);
        fflush(stdout);

        for (i = 1; i < npolys; i++)
        {
            if (!fmpz_mpoly_get_fmpz_poly(res + 0, mres + i, 0, ctx) ||
                !fmpz_poly_equal(res + 0, res + i))
            {
                flint_printf("oops");
                flint_abort();
            }
        }
    }
}

int main(int argc, char *argv[])
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("\n");

    {
        slong npolys = 10000;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_struct * mpolys, * mres;
        fmpz_poly_struct * polys, * res;
        slong len;
        ulong exp_bound;
        flint_bitcnt_t coeff_bits;

        fmpz_mpoly_ctx_init(ctx, 1, ORD_LEX);
        mpolys = FLINT_ARRAY_ALLOC(npolys, fmpz_mpoly_struct);
        polys = FLINT_ARRAY_ALLOC(npolys, fmpz_poly_struct);
        mres = FLINT_ARRAY_ALLOC(npolys, fmpz_mpoly_struct);
        res = FLINT_ARRAY_ALLOC(npolys, fmpz_poly_struct);

        for (i = 0; i < npolys; i++)
        {
            fmpz_mpoly_init(mpolys + i, ctx);
            fmpz_poly_init(polys + i);
            fmpz_mpoly_init(mres + i, ctx);
            fmpz_poly_init(res + i);
        }

        flint_printf("****** dense mul, small coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 10);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 10);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** dense mul, large coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 10);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 200);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** large dense mul, small coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 30);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 10);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** large dense mul, large coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 30);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 200);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** sparse mul, small coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 10);
            exp_bound = n_randint(state, 10*len) + 10;
            coeff_bits = 10 + n_randint(state, 10);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** sparse mul, large coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 10);
            exp_bound = n_randint(state, 10*len) + 10;
            coeff_bits = 10 + n_randint(state, 100);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** very sparse mul, small coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 10);
            exp_bound = n_randint(state, 20*len) + 10;
            coeff_bits = 10 + n_randint(state, 10);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** very sparse mul, large coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 10);
            exp_bound = n_randint(state, 20*len) + 10;
            coeff_bits = 10 + n_randint(state, 100);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, fmpz_mpoly_mul, fmpz_poly_mul);

        flint_printf("****** gcd, small coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 5);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 10);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
        }
        for (i = 0; i + 2 <= npolys; i += 2)
        {
            len = 1 + n_randint(state, 5);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 10);
            fmpz_mpoly_randtest_bound(mres + 0, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_mul(mpolys + i, mpolys + i, mres + 0, ctx);
            fmpz_mpoly_mul(mpolys + i + 1, mpolys + i + 1, mres + 0, ctx);
        }
        for (i = 0; i < npolys; i++)
        {
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, (void (*)(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t)) fmpz_mpoly_gcd, fmpz_poly_gcd);

        flint_printf("****** gcd, large coeffs ********\n");
        for (i = 0; i < npolys; i++)
        {
            len = 1 + n_randint(state, 5);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 200);
            fmpz_mpoly_randtest_bound(mpolys + i, state, len, coeff_bits, exp_bound, ctx);
        }
        for (i = 0; i + 2 <= npolys; i += 2)
        {
            len = 1 + n_randint(state, 5);
            exp_bound = len + 1;
            coeff_bits = 10 + n_randint(state, 200);
            fmpz_mpoly_randtest_bound(mres + 0, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_mul(mpolys + i, mpolys + i, mres + 0, ctx);
            fmpz_mpoly_mul(mpolys + i + 1, mpolys + i + 1, mres + 0, ctx);
        }
        for (i = 0; i < npolys; i++)
        {
            fmpz_mpoly_get_fmpz_poly(polys + i, mpolys + i, 0, ctx);
        }
        run_cmp(mpolys, polys, mres, res, npolys, ctx, (void (*)(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t)) fmpz_mpoly_gcd, fmpz_poly_gcd);

        for (i = 0; i < npolys; i++)
        {
            fmpz_mpoly_clear(mpolys + i, ctx);
            fmpz_poly_clear(polys + i);
            fmpz_mpoly_clear(mres + i, ctx);
            fmpz_poly_clear(res + i);
        }
        flint_free(mpolys);
        flint_free(polys);
        flint_free(mres);
        flint_free(res);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
