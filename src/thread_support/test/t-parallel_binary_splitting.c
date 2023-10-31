/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "fmpz.h"

typedef struct
{
    fmpz r;
}
product_res_t;

typedef struct
{
    mp_srcptr factors;
    int left_inplace;
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
    slong i;

    fmpz_one(&res->r);

    for (i = a; i < b; i++)
        fmpz_mul_ui(&res->r, &res->r, args->factors[i]);
}

static void
bsplit_product(fmpz_t r, mp_srcptr factors, slong len, slong thread_limit, int flags)
{
    product_res_t res;
    product_args_t args;

    res.r = *r;

    args.factors = factors;
    args.left_inplace = (flags & FLINT_PARALLEL_BSPLIT_LEFT_INPLACE) ? 1 : 0;

    flint_parallel_binary_splitting(&res,
        (bsplit_basecase_func_t) product_basecase,
        (bsplit_merge_func_t) product_combine,
        sizeof(product_res_t),
        (bsplit_init_func_t) product_init,
        (bsplit_clear_func_t) product_clear,
        &args, 0, len, 4, thread_limit, flags);

    *r = res.r;
}

TEST_FUNCTION_START(thread_support_parallel_binary_splitting, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_t r, s;
        mp_ptr factors;
        slong i, n;
        int flags;

        n = n_randint(state, 100);

        flint_set_num_threads(n_randint(state, 10) + 1);

        factors = flint_malloc(n * sizeof(mp_limb_t));

        fmpz_init(r);
        fmpz_init(s);

        for (i = 0; i < n; i++)
            factors[i] = n_randint(state, 300);

        flags = 0;
        if (n_randint(state, 2))
            flags = FLINT_PARALLEL_BSPLIT_LEFT_INPLACE;

        bsplit_product(r, factors, n, n_randint(state, 5), flags);

        fmpz_one(s);
        for (i = 0; i < n; i++)
            fmpz_mul_ui(s, s, factors[i]);

        if (!fmpz_equal(r, s))
        {
            flint_printf("FAIL\n");
            flint_printf("num_threads = %wd, i = %wd/%wd\n", flint_get_num_threads(), i, n);
            flint_abort();
        }

        flint_free(factors);
        fmpz_clear(r);
        fmpz_clear(s);
    }

    TEST_FUNCTION_END(state);
}
