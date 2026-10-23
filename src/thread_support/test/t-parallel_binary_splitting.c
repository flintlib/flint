/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "fmpz.h"

typedef struct
{
    fmpz r;
}
product_res_t;

typedef struct
{
    nn_srcptr factors;
    int left_inplace;
    slong max_workers;
    int too_many_workers;
}
product_args_t;

static void
product_init(product_res_t * x, product_args_t * args)
{
    fmpz_init(&x->r);
}

static void
product_clear(product_res_t * x, product_args_t * args)
{
    fmpz_clear(&x->r);
}

static void
product_combine(product_res_t * res, product_res_t * left, product_res_t * right, product_args_t * args)
{
    if (((res == left) != args->left_inplace) || res == right)
    {
        flint_abort();
    }

    fmpz_mul(&res->r, &left->r, &right->r);
}

static void
product_basecase(product_res_t * res, slong a, slong b, product_args_t * args)
{
    slong i, held;

    /* Workers of the global pool held by this computation (the test does
       not use the pool otherwise). This number only changes when a split
       requests or gives back a worker, so the check is deterministic:
       each split holds its worker until its whole subtree is done. */
    held = thread_pool_get_size(global_thread_pool)
               - (flint_get_num_available_threads() - 1);
    if (held > args->max_workers)
        args->too_many_workers = 1;

    fmpz_one(&res->r);

    for (i = a; i < b; i++)
        fmpz_mul_ui(&res->r, &res->r, args->factors[i]);
}

/* returns 1 if the computation used more workers than allowed */
static int
bsplit_product(fmpz_t r, nn_srcptr factors, slong len, slong thread_limit, int flags)
{
    product_res_t res;
    product_args_t args;

    res.r = *r;

    args.factors = factors;
    args.left_inplace = (flags & FLINT_PARALLEL_BSPLIT_LEFT_INPLACE) ? 1 : 0;
    args.max_workers = flint_get_num_threads() - 1;
    if (thread_limit > 0)
        args.max_workers = FLINT_MIN(args.max_workers, thread_limit - 1);
    args.too_many_workers = 0;

    flint_parallel_binary_splitting(&res,
        (bsplit_basecase_func_t) product_basecase,
        (bsplit_merge_func_t) product_combine,
        sizeof(product_res_t),
        (bsplit_init_func_t) product_init,
        (bsplit_clear_func_t) product_clear,
        &args, 0, len, 4, thread_limit, flags);

    *r = res.r;

    return args.too_many_workers;
}

TEST_FUNCTION_START(thread_support_parallel_binary_splitting, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_t r, s;
        nn_ptr factors;
        slong i, n, thread_limit;
        int flags, too_many_workers;

        n = n_randint(state, 100);

        flint_set_num_threads(n_randint(state, 10) + 1);

        /* the case of flint_parallel_binary_splitting ignoring thread_limit */
        if (iter == 0)
        {
            n = 64;
            flint_set_num_threads(8);
        }

        factors = flint_malloc(n * sizeof(ulong));

        fmpz_init(r);
        fmpz_init(s);

        for (i = 0; i < n; i++)
            factors[i] = n_randint(state, 300);

        flags = 0;
        if (n_randint(state, 2))
            flags = FLINT_PARALLEL_BSPLIT_LEFT_INPLACE;

        thread_limit = (iter == 0) ? 3 : n_randint(state, 5);
        too_many_workers = bsplit_product(r, factors, n, thread_limit, flags);

        fmpz_one(s);
        for (i = 0; i < n; i++)
            fmpz_mul_ui(s, s, factors[i]);

        if (!fmpz_equal(r, s))
            TEST_FUNCTION_FAIL(
                    "num_threads = %wd, i = %wd/%wd\n",
                    flint_get_num_threads(), i, n);

        if (too_many_workers)
            TEST_FUNCTION_FAIL(
                    "too many workers: num_threads = %wd, thread_limit = %wd, n = %wd\n",
                    flint_get_num_threads(), thread_limit, n);

        flint_free(factors);
        fmpz_clear(r);
        fmpz_clear(s);
    }

    TEST_FUNCTION_END(state);
}
