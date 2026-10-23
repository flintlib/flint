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
#include "mpn_extras.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* The shared tables under concurrent growth: several threads evaluate
   the table-driven functions at random precisions from cold caches,
   each reading whatever version of the tables is current and racing
   the others to grow them, and the results are compared with a serial
   rerun of the same evaluations against the final tables.  Every
   evaluation reads the top limbs or digits of its entries, which are
   the same floors in every version, so the results must agree
   bitwise. */

#define NTHREADS 8
#define NEVALS 40
#define NMAX 160
#define OUTLEN (2 * (NMAX + 2))

typedef struct
{
    ulong * out;        /* NTHREADS * NEVALS * OUTLEN limbs */
    ulong seed;
}
work_struct;

/* the evaluations of worker i, in order */
static void
evaluate(ulong * out, slong i, ulong seed)
{
    flint_rand_t state;
    slong iter, n, k;
    ulong x[NMAX], * y, * y2;

    flint_rand_init(state);
    flint_rand_set_seed(state, seed + 1000 * i, 0);

    for (iter = 0; iter < NEVALS; iter++, out += OUTLEN)
    {
        int which = n_randint(state, 9);
        int r;

        flint_mpn_zero(out, OUTLEN);
        y = out;
        y2 = out + NMAX + 2;

        n = 1 + n_randint(state, NMAX);
        /* the diophantine reductions and the constants from 2 limbs,
           the 32-bit bitwise paths too */
        if (which >= 5 || FLINT_BITS == 32)
            n = FLINT_MAX(n, 2);

        for (k = 0; k < n; k++)
            x[k] = n_randtest(state);
        /* arguments in [0, 1), away from 1 for the trigonometric
           functions' halving to stay in range */
        x[n - 1] >>= 1;

        r = (n_randint(state, 3) == 0)
            ? 32 + (int) n_randint(state, FLINT_BITS * n - 48 + 1) : 0;
        if (n == 1 && r != 0)
            r = 32;

        switch (which)
        {
            case 0: _mp_real_exp_bitwise_rs(y, NULL, x, n, r); break;
            case 1: _mp_real_log1p_bitwise_rs(y, NULL, x, n, r); break;
            case 2: _mp_real_atan_bitwise_rs(y, NULL, x, n, r); break;
            case 3: _mp_real_sin_cos_bitwise_rs(y, y2, NULL, x, n, r); break;
            case 4: _mp_real_tan_bitwise_rs(y, NULL, x, n, r); break;
            case 5: _mp_real_exp_diophantine(y, NULL, x, n); break;
            case 6: _mp_real_sin_cos_diophantine(y, y2, NULL, x, n); break;
            case 7: _mp_real_const_pi4(y, NULL, n, 1); break;
            default: _mp_real_const_log2(y, NULL, n, 1); break;
        }
    }

    flint_rand_clear(state);
}

static void
worker(slong i, void * args)
{
    work_struct * W = args;

    evaluate(W->out + i * NEVALS * OUTLEN, i, W->seed);
}

static void
clear_all(void)
{
#if FLINT_BITS == 64
    _mp_real_bitwise_vec_clear();
#endif
    _mp_real_exp_logs_clear();
    _mp_real_atans_clear();
    _mp_real_log_primes_clear();
    _mp_real_atan_gauss_clear();
    _mp_real_const_clear_cache();
}

TEST_FUNCTION_START(mp_real_shared_tables, state)
{
    slong iter, i;
    int old_threads = flint_get_num_threads();

    for (iter = 0; iter < 1 + flint_test_multiplier() / 2; iter++)
    {
        work_struct W;
        ulong * ref;

        W.out = flint_malloc(NTHREADS * NEVALS * OUTLEN * sizeof(ulong));
        W.seed = n_randtest(state);
        ref = flint_malloc(NEVALS * OUTLEN * sizeof(ulong));

        clear_all();
        flint_set_num_threads(NTHREADS);
        flint_parallel_do(worker, &W, NTHREADS, NTHREADS, FLINT_PARALLEL_UNIFORM);
        flint_set_num_threads(old_threads);

        for (i = 0; i < NTHREADS; i++)
        {
            evaluate(ref, i, W.seed);
            if (mpn_cmp(ref, W.out + i * NEVALS * OUTLEN, NEVALS * OUTLEN) != 0)
                TEST_FUNCTION_FAIL("worker %wd: results differ from the serial rerun\n", i);
        }

        flint_free(W.out);
        flint_free(ref);
    }

    clear_all();

    TEST_FUNCTION_END(state);
}
