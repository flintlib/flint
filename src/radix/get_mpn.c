/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"
#include "thread_pool.h"
#include "thread_support.h"

static void
radix_powers_init_ui(radix_powers_t powers, ulong b, ulong max_exp)
{
    slong k, len, n, val, bound_limbs;
    nn_ptr bpow, bbuf;

    FLINT_ASSERT(b != 0);
    FLINT_ASSERT(max_exp >= 1);

    len = 1;
    powers->exps[0] = max_exp;
    bound_limbs = max_exp + 1;

    while (powers->exps[len - 1] > 1)
    {
        powers->exps[len] = (powers->exps[len - 1] + 1) / 2;
        bound_limbs += powers->exps[len] + 1;
        len++;
    }
    powers->len = len;

    powers->buf = bbuf = flint_malloc(bound_limbs * sizeof(ulong));

    /* pows[len - 1] = b ^ 1 */
    powers->pows[len - 1] = bbuf;
    powers->pows[len - 1][0] = b;
    powers->sizes[len - 1] = 1;
    powers->val_limbs[len - 1] = 0;
    bbuf += 1;

    for (k = len - 2; k >= 0; k--)
    {
        n = powers->sizes[k + 1];
        val = powers->val_limbs[k + 1];
        powers->pows[k] = bbuf;
        bpow = powers->pows[k];

        flint_mpn_sqr(bpow + 2 * val, powers->pows[k + 1] + val, n - val);
        flint_mpn_zero(bpow, 2 * val);
        n = 2 * n;
        n -= (bpow[n - 1] == 0);

        if (powers->exps[k] != 2 * powers->exps[k + 1])
        {
            mpn_divexact_1(bpow, bpow, n, b);
            n -= (bpow[n - 1] == 0);

            val = 0;
            while (bpow[val] == 0)
                val++;
        }
        else
        {
            val = 2 * val;
            val += (bpow[val] == 0);
        }

        bbuf += n;

        powers->sizes[k] = n;
        powers->val_limbs[k] = val;

        FLINT_ASSERT(val == 0 || mpn_zero_p(bpow, val));
    }
}

void
radix_powers_clear(radix_powers_t powers)
{
    flint_free(powers->buf);
}

slong
radix_get_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    ulong B = LIMB_RADIX(radix), cy, lo;
    slong i, len;

    while (an > 0 && a[an - 1] == 0)
        an--;

    if (an == 0)
        return 0;

    if (an == 1)
    {
        res[0] = a[0];
        return 1;
    }

    umul_ppmm(cy, lo, a[an - 1], B);
    add_ssaaaa(cy, lo, cy, lo, 0, a[an - 2]);
    res[0] = lo;
    len = 1;
    if (cy != 0)
    {
        res[1] = cy;
        len = 2;
    }

    for (i = an - 3; i >= 0; i--)
    {
        cy = mpn_mul_1(res, res, len, B);
        cy += mpn_add_1(res, res, len, a[i]);
        if (cy != 0)
        {
            res[len] = cy;
            len++;
        }
    }

    return len;
}

static slong
_radix_get_mpn_recursive(nn_ptr res, nn_srcptr a, slong an, const radix_powers_t powers, slong depth, const radix_t radix);

typedef struct
{
    nn_ptr res;
    nn_srcptr a;
    slong an;
    const radix_powers_struct * powers;
    slong depth;
    const radix_struct * radix;
    slong out_len;
}
worker_args_struct;

static void
worker(void * arg)
{
    worker_args_struct * X = (worker_args_struct * ) arg;
    X->out_len = _radix_get_mpn_recursive(X->res, X->a, X->an, X->powers, X->depth, X->radix);
}

static slong
_radix_get_mpn_recursive(nn_ptr res, nn_srcptr a, slong an, const radix_powers_t powers, slong depth, const radix_t radix)
{
    slong an_lo, an_hi;
    nn_ptr lo, hi;
    slong len_lo, len_hi, len_pow, len_res;
    ulong cy;
    TMP_INIT;

    while (an > 0 && a[an - 1] == 0)
        an--;

    if (an < 50)
        return radix_get_mpn_basecase(res, a, an, radix);

    while (powers->exps[depth] > (an + 1) / 2)
        depth++;

    an_lo = powers->exps[depth];
    an_hi = an - an_lo;

    /* todo: can we tighten these allocations? */
    TMP_START;
    lo = TMP_ALLOC((an + 2) * sizeof(ulong));
    hi = lo + (an_lo + 1);

    slong nworkers, nthreads, nworkers_save;
    int want_workers;
    thread_pool_handle * threads;

    nworkers = 0;
    if (an > 500)
    {
        nthreads = flint_get_num_threads();
        /* Prefer to let the multithreaded
           multiplication do its things near the root. */
        want_workers = nthreads >= 2 && (an_lo <= 100000000 || depth >= 2);
        nworkers = flint_request_threads(&threads, want_workers ? 2 : 1);
    }

    if (nworkers == 1)
    {
        worker_args_struct work_lo = { lo, a, an_lo, powers, depth, radix, 0 };
        worker_args_struct work_hi = { hi, a + an_lo, an_hi, powers, depth, radix, 0 };

        nworkers_save = flint_set_num_workers(nthreads - nthreads / 2 - 1);
        thread_pool_wake(global_thread_pool, threads[0], nthreads / 2 - 1, worker, &work_lo);
        worker(&work_hi);

        flint_reset_num_workers(nworkers_save);
        thread_pool_wait(global_thread_pool, threads[0]);
        flint_give_back_threads(threads, nworkers);

        len_lo = work_lo.out_len;
        len_hi = work_hi.out_len;
    }
    else
    {
        len_lo = _radix_get_mpn_recursive(lo, a, an_lo, powers, depth, radix);
        len_hi = _radix_get_mpn_recursive(hi, a + an_lo, an_hi, powers, depth, radix);
    }

    if (len_hi == 0)
    {
        flint_mpn_copyi(res, lo, len_lo);
        len_res = len_lo;
    }
    else
    {
        slong val = powers->val_limbs[depth];
        len_pow = powers->sizes[depth] - val;

        if (len_hi >= len_pow)
            flint_mpn_mul(res + val, hi, len_hi, powers->pows[depth] + val, len_pow);
        else
            flint_mpn_mul(res + val, powers->pows[depth] + val, len_pow, hi, len_hi);

        len_res = len_pow + len_hi + val;
        len_res -= (res[len_res - 1] == 0);

        FLINT_ASSERT(len_res >= len_lo);

        if (len_lo <= val)
        {
            flint_mpn_copyi(res, lo, len_lo);
            flint_mpn_zero(res + len_lo, val - len_lo);
        }
        else
        {
            flint_mpn_copyi(res, lo, val);
            cy = mpn_add(res + val, res + val, len_res - val, lo + val, len_lo - val);
            if (cy != 0)
            {
                res[len_res] = cy;
                len_res++;
            }
        }
    }

    TMP_END;

    return len_res;
}

slong
radix_get_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    slong rn;
    radix_powers_t powers;
    radix_powers_init_ui(powers, LIMB_RADIX(radix), FLINT_MAX((an + 1) / 2, 1));
    rn = _radix_get_mpn_recursive(res, a, an, powers, 0, radix);
    radix_powers_clear(powers);
    return rn;
}

slong
radix_get_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    FLINT_ASSERT(res != a);

    while (an > 0 && a[an - 1] == 0)
        an--;

    if (an < 50)
        return radix_get_mpn_basecase(res, a, an, radix);
    else
        return radix_get_mpn_divconquer(res, a, an, radix);
}

