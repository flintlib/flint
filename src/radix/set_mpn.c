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
radix_powers_init_ui_radix(radix_powers_t powers, ulong b, ulong max_exp, const radix_t radix)
{
    slong k, len, n, bound_limbs;
    nn_ptr bpow, bbuf;

    FLINT_ASSERT(b != 0);
    FLINT_ASSERT(max_exp >= 1);
    FLINT_ASSERT(b < LIMB_RADIX(radix));

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
        powers->pows[k] = bbuf;
        bpow = powers->pows[k];

        radix_sqr(bpow, powers->pows[k + 1], n, radix);
        n = 2 * n;
        n -= (bpow[n - 1] == 0);

        if (powers->exps[k] != 2 * powers->exps[k + 1])
        {
            radix_divexact_1(bpow, bpow, n, b, radix);
            n -= (bpow[n - 1] == 0);
        }

        bbuf += n;
        powers->sizes[k] = n;
        powers->val_limbs[k] = 0;
    }
}

/* Extract bits [start, stop), writing ceil((stop-start)/FLINT_BITS) limbs. */
static void
flint_mpn_slice_bits(nn_ptr res, nn_srcptr x, mp_bitcnt_t start, mp_bitcnt_t stop)
{
    FLINT_ASSERT(start >= 0);
    FLINT_ASSERT(start < stop);

    slong res_limbs = (stop - start + FLINT_BITS - 1) / FLINT_BITS;
    slong start_limb = start / FLINT_BITS;
    slong stop_limb = stop / FLINT_BITS;
    slong rem;

    rem = start % FLINT_BITS;

    if (rem == 0)
    {
        mpn_copyi(res, x + start_limb, res_limbs);
    }
    else
    {
        mpn_rshift(res, x + start_limb, res_limbs, rem);

        if (stop_limb == start_limb + res_limbs && stop_limb * FLINT_BITS < stop)
            res[res_limbs - 1] |= (x[stop_limb] << (FLINT_BITS - rem));
    }

    rem = (stop - start) % FLINT_BITS;
    if (rem != 0)
        res[res_limbs - 1] &= ((UWORD(1) << rem) - 1);
}

slong
radix_set_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    ulong B = LIMB_RADIX(radix);
    slong len;
    nn_ptr tmp;
    TMP_INIT;

    while (an > 0 && a[an - 1] == 0)
        an--;

    if (an == 0)
        return 0;

    if (an == 1 && a[0] < B)
    {
        res[0] = a[0];
        return 1;
    }

    TMP_START;
    tmp = TMP_ALLOC(an * sizeof(ulong));

    flint_mpn_copyi(tmp, a, an);

    len = 0;

    while (an > 12)
    {
        res[len] = mpn_divrem_1(tmp, 0, tmp, an, B);
        len++;
        an -= (tmp[an - 1] == 0);
    }

    while (an != 0)
    {
        res[len] = flint_mpn_divrem_1_preinv(tmp, tmp, an, radix->B.n, radix->B.ninv, radix->B.norm);
        len++;
        an -= (tmp[an - 1] == 0);
    }

    TMP_END;

    return len;
}

static slong
radix_set_mpn_basecase1(nn_ptr res, nn_srcptr a, slong an_words, slong abits, slong astart, slong astop, const radix_t radix)
{
    slong tbits, tlimbs, len;
    nn_ptr tmp;
    TMP_INIT;

    astart *= abits;
    astop *= abits;

    astop = FLINT_MIN(astop, an_words * FLINT_BITS);

    tbits = astop - astart;
    tlimbs = (tbits + FLINT_BITS - 1) / FLINT_BITS;

    TMP_START;
    tmp = TMP_ALLOC(tlimbs * sizeof(ulong));
    flint_mpn_slice_bits(tmp, a, astart, astop);
    len = radix_set_mpn_basecase(res, tmp, tlimbs, radix);

    TMP_END;

    return len;
}

static slong
_radix_set_mpn_recursive(nn_ptr res, nn_srcptr a, slong an_words, slong abits, slong astart, slong astop, const radix_powers_t powers, slong depth, const radix_t radix);

typedef struct
{
    nn_ptr res;
    nn_srcptr a;
    slong an_words;
    slong abits;
    slong astart;
    slong astop;
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
    X->out_len = _radix_set_mpn_recursive(X->res, X->a, X->an_words, X->abits, X->astart, X->astop, X->powers, X->depth, X->radix);
}

static slong
_radix_set_mpn_recursive(nn_ptr res, nn_srcptr a, slong an_words, slong abits, slong astart, slong astop, const radix_powers_t powers, slong depth, const radix_t radix)
{
    slong an, an_lo;
    nn_ptr lo, hi;
    slong len_lo, len_hi, len_pow, len_res;
    ulong cy;
    TMP_INIT;

    an = astop - astart;

    if (an < 25)
    {
        len_lo = radix_set_mpn_basecase1(res, a, an_words, abits, astart, astop, radix);
        return len_lo;
    }

    while (powers->exps[depth] > (an + 1) / 2)
        depth++;

    an_lo = powers->exps[depth];

    /* todo: can we tighten these allocations? */
    TMP_START;
    lo = TMP_ALLOC((an + 2) * sizeof(ulong));
    hi = lo + an_lo + 1;

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
        worker_args_struct work_lo = { lo, a, an_words, abits, astart, astart + an_lo, powers, depth, radix, 0 };
        worker_args_struct work_hi = { hi, a, an_words, abits, astart + an_lo, astop, powers, depth, radix, 0 };

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
        len_lo = _radix_set_mpn_recursive(lo, a, an_words, abits, astart, astart + an_lo, powers, depth, radix);
        len_hi = _radix_set_mpn_recursive(hi, a, an_words, abits, astart + an_lo, astop, powers, depth, radix);
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
            radix_mul(res + val, hi, len_hi, powers->pows[depth] + val, len_pow, radix);
        else
            radix_mul(res + val, powers->pows[depth] + val, len_pow, hi, len_hi, radix);

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
            cy = radix_add(res + val, res + val, len_res - val, lo + val, len_lo - val, radix);
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
radix_set_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    FLINT_ASSERT(res != a);

    slong abits, abits_limbs;
    slong rn;
    radix_powers_t powers;
    ulong B = LIMB_RADIX(radix);

    /* As a special case, the divide and conquer code will not work as intended
       for radix B = 2. Should put in special-purpose code for all binary
       radices; this is just a quick workaround. */
    if (B == 2)
    {
        slong i, j;
        while (an > 0 && a[an - 1] == 0)
            an--;

        for (i = 0; i < an; i++)
        {
            ulong c = a[i];

            for (j = 0; j < FLINT_BITS; j++)
                res[i * FLINT_BITS + j] = ((c >> j) & UWORD(1));
        }

        rn = an * FLINT_BITS;
        while (rn > 0 && res[rn - 1] == 0)
            rn--;
        return rn;
    }

    abits = FLINT_BIT_COUNT(B - 1) - 1;
    abits_limbs = (an * FLINT_BITS + abits - 1) / abits;

    slong e = FLINT_MAX((abits_limbs + 1) / 2, 1);

    radix_powers_init_ui_radix(powers, UWORD(1) << abits, e, radix);
    rn = _radix_set_mpn_recursive(res, a, an, abits, 0, abits_limbs, powers, 0, radix);
    radix_powers_clear(powers);

    return rn;
}

slong
radix_set_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    FLINT_ASSERT(res != a);

    if (an < 40)
        return radix_set_mpn_basecase(res, a, an, radix);
    else
        return radix_set_mpn_divconquer(res, a, an, radix);
}

slong
radix_set_mpn_need_alloc(slong n, const radix_t radix)
{
    if (n == 0)
        return 0;

    slong Bbits = NMOD_BITS(radix->B) - 1;
    slong nlimbs = (n * FLINT_BITS + Bbits - 1) / Bbits;
    return nlimbs + 1;
}

