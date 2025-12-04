/*
    Copyright (C) 2015 Kushagra Singh
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"

static ulong
n_sqr_and_add_a_fast(ulong y, ulong a, nmod_redc_ctx_t ctx)
{
    y = nmod_redc_fast_mul(y, y, ctx);
    /* Has a negligible chance of not being reduced, which is fine. */
    y += a;
    return y;
}

static ulong
n_sqr_and_add_a(ulong y, ulong a, nmod_redc_ctx_t ctx)
{
    y = nmod_redc_mul(y, y, ctx);
    /* Has a negligible chance of not being reduced, which is fine. */
    y += a;
    return y;
}

static ulong
n_absdiff(ulong x, ulong y)
{
    if (x > y)
        return x - y;
    else
        return y - x;
}

int
n_factor_pollard_brent_single(ulong *factor, ulong n,
                              ulong ai, ulong xi, ulong max_iters)
{
    ulong iter, i, k, j, minval, m, x, y, a, q, ys;
    int ret;
    int fast;

    if (n < 4)
        return 0;

    nmod_redc_ctx_t ctx;
    nmod_redc_ctx_init_ui(ctx, n);
    fast = nmod_redc_can_use_fast(ctx);

    q = 1;
    (*factor) = 1;

    y = xi;
    a = ai;

    m = 100;
    iter = 1;
    do {
        x = y;
        k = 0;

        if (fast)
        {
            for (i = 0; i < iter; i++)
                y = n_sqr_and_add_a_fast(y, a, ctx);
        }
        else
        {
            for (i = 0; i < iter; i++)
                y = n_sqr_and_add_a(y, a, ctx);
        }

        do {
            minval = iter - k;
            if (m < minval)
                minval = m;

            ys = y;

            if (fast)
            {
                for (i = 0; i < minval; i++)
                {
                    y = n_sqr_and_add_a_fast(y, a, ctx);
                    q = nmod_redc_fast_mul(q, n_absdiff(x, y), ctx);
                }
            }
            else
            {
                for (i = 0; i < minval; i++)
                {
                    y = n_sqr_and_add_a(y, a, ctx);
                    q = nmod_redc_mul(q, n_absdiff(x, y), ctx);
                }
            }

            if (q == 0)
                (*factor) = n;
            else
                (*factor) = n_gcd(q, n);

            k += m;
            j = ((*factor) == 1);
        } while ((k < iter) && (j));

        if (iter > max_iters)
            break;
        iter *= 2;
    }  while (j);

    if ((*factor) == n)
    {
        if (fast)
        {
            do {
                ys = n_sqr_and_add_a_fast(ys, a, ctx);
                if (q == 0)
                    (*factor) = n;
                else
                    (*factor) = n_gcd(q, n);
                (*factor) = n_gcd(n_absdiff(x, ys), n);
            } while ((*factor) == 1);   /* gcd == 1 */
        }
        else
        {
            do {
                ys = n_sqr_and_add_a(ys, a, ctx);
                if (q == 0)
                    (*factor) = n;
                else
                    (*factor) = n_gcd(q, n);
                (*factor) = n_gcd(n_absdiff(x, ys), n);
            } while ((*factor) == 1);   /* gcd == 1 */
        }
    }

    ret = 1;

    if ((*factor) == 1) /* gcd == 1 */
        ret = 0;
    else if ((*factor) == n) /* gcd == n*/
        ret = 0;

    return ret;
}

int
n_factor_pollard_brent(ulong *factor, flint_rand_t state, ulong n,
                        ulong max_tries, ulong max_iters)
{
    ulong a, x;
    int ret = 0;

    while (max_tries--)
    {
        /* Keep a small so that additions without reduction have a high
           chance of remaining valid. */
        a = n_randint(state, FLINT_MIN(n - 4, 1024)) + 1;  /* 1 <= a <= n - 3 */
        x = n_randint(state, n - 2) + 1;  /* 1 <= a <= n - 3 */

        ret = n_factor_pollard_brent_single(factor, n, a, x, max_iters);

        if (ret == 1)
            return 1;

        a += 1;
    }

    return ret;
}
