/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <ctype.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "thread_support.h"
#include "thread_pool.h"
#include "profiler.h"

#define BASECASE_CUTOFF 24000

/* lower bounds */
#if FLINT_BITS == 64
#define DIGITS_PER_LIMB 19
#else
#define DIGITS_PER_LIMB 9
#endif

/* todo: the binary splitting code is almost the same as that for
   get_str and could be reused */
typedef struct
{
    fmpz * res;
    const char * s;
    slong slen;
    const slong * exps;
    slong cur_depth;
    slong depth;
    const fmpz * pows;
    const fmpz_preinvn_struct * preinv;
}
worker_args_struct;

static void
_fmpz_get_str_recursive(fmpz_t res, const char * s, slong slen, const slong * exps, slong cur_depth, slong depth, const fmpz * pows);

static void
worker(void * arg)
{
    worker_args_struct * X = (worker_args_struct * ) arg;
    _fmpz_get_str_recursive(X->res, X->s, X->slen, X->exps, X->cur_depth, X->depth, X->pows);
}

static void
_fmpz_set_str_basecase(fmpz_t res, const char * s, slong slen)
{
    mp_ptr tmp;
    unsigned char * stmp;
    mp_size_t n;
    slong i;
    TMP_INIT;

    TMP_START;

    stmp = TMP_ALLOC(sizeof(char) * slen);
    tmp = TMP_ALLOC(sizeof(mp_limb_t) * (slen / DIGITS_PER_LIMB + 2));

    for (i = 0; i < slen; i++)
        stmp[i] = s[i] - '0';

    n = mpn_set_str(tmp, stmp, slen, 10);

    if (n == 0)
        fmpz_zero(res);
    else
        fmpz_set_ui_array(res, tmp, n);

    TMP_END;
}

static void
_fmpz_get_str_recursive(fmpz_t res, const char * s, slong slen, const slong * exps, slong cur_depth, slong depth, const fmpz * pows)
{
    if (cur_depth >= depth || slen < BASECASE_CUTOFF)
    {
        _fmpz_set_str_basecase(res, s, slen);
    }
    else
    {
        fmpz_t q, r;
        slong num_right = exps[cur_depth];
        slong nworkers, nthreads, nworkers_save;
        int want_workers;
        thread_pool_handle * threads;
        worker_args_struct high_digits[1], low_digits[1];

        fmpz_init(q);
        fmpz_init(r);

        /*
        _fmpz_get_str_recursive(r, s + slen - num_right, num_right, exps, cur_depth + 1, depth, pows);
        _fmpz_get_str_recursive(q, s, slen - num_right, exps, cur_depth + 1, depth, pows);
        fmpz_mul(res, q, pows + cur_depth);
        fmpz_mul_2exp(res, res, exps[cur_depth]);
        fmpz_add(res, res, r);
        */

        low_digits->res = r;
        low_digits->s = s + slen - num_right;
        low_digits->slen = num_right;
        low_digits->exps = exps;
        low_digits->cur_depth = cur_depth + 1;
        low_digits->depth = depth;
        low_digits->pows = pows;

        high_digits->res = q;
        high_digits->s = s;
        high_digits->slen = slen - num_right;
        high_digits->exps = exps;
        high_digits->cur_depth = cur_depth + 1;
        high_digits->depth = depth;
        high_digits->pows = pows;

        nthreads = flint_get_num_threads();

        /* Prefer to let the multithreaded
           multiplication do its things near the root. */
        want_workers = nthreads >= 2 && (num_right <= 100000000 || cur_depth >= 2);
        nworkers = flint_request_threads(&threads, want_workers ? 2 : 1);

        if (nworkers == 1)
        {
            nworkers_save = flint_set_num_workers(nthreads - nthreads / 2 - 1);
            thread_pool_wake(global_thread_pool, threads[0], nthreads / 2 - 1, worker, low_digits);
            worker(high_digits);
            flint_reset_num_workers(nworkers_save);
            thread_pool_wait(global_thread_pool, threads[0]);
        }
        else
        {
            worker(low_digits);
            worker(high_digits);
        }

        flint_give_back_threads(threads, nworkers);

        fmpz_mul(res, q, pows + cur_depth);
        fmpz_mul_2exp(res, res, exps[cur_depth]);
        fmpz_add(res, res, r);

        fmpz_clear(q);
        fmpz_clear(r);
    }
}

void
fmpz_set_str_bsplit_threaded(fmpz_t res, const char * s, slong slen)
{
    slong k, depth;
    slong exps[FLINT_BITS];
    fmpz * pows;

    exps[0] = (slen + 1) / 2;
    depth = 1;

    while (exps[depth - 1] > BASECASE_CUTOFF / DIGITS_PER_LIMB)
    {
        exps[depth] = (exps[depth - 1] + 1) / 2;
        depth++;
    }

    pows = _fmpz_vec_init(depth);

    fmpz_ui_pow_ui(pows + depth - 1, 5, exps[depth - 1]);
    for (k = depth - 2; k >= 0; k--)
    {
        fmpz_mul(pows + k, pows + k + 1, pows + k + 1);
        if (exps[k] != 2 * exps[k + 1])
            fmpz_divexact_ui(pows + k, pows + k, 5);
    }

    _fmpz_get_str_recursive(res, s, slen, exps, 0, depth, pows);
    _fmpz_vec_clear(pows, depth);
}

static int
fmpz_set_str_fallback(fmpz_t res, const char * str, int b, int neg)
{
    int err;
    __mpz_struct * z = _fmpz_promote(res);
    err = mpz_set_str(z, str, b);
    if (neg)
        mpz_neg(z, z);
    _fmpz_demote_val(res);
    return err;
}

int
fmpz_set_str(fmpz_t res, const char * str, int base)
{
    slong slen, i;
    int neg = 0;

    /* Let GMP handle unusual bases. */
    if (base != 10)
        return fmpz_set_str_fallback(res, str, base, 0);

    /* Allow leading whitespace. */
    while (isspace(str[0]))
        str++;

    if (str[0] == '-')
    {
        str++;
        neg = 1;
    }

    slen = strlen(str);

    /* Allow trailing whitespace. */
    while (slen > 0 && isspace(str[slen - 1]))
        slen--;

    if (slen == 0)
        return -1;

    for (i = 0; i < slen; i++)
    {
        /* The string is either invalid or has interior whitespace,
           which GMP allows. Either way, let GMP handle it. */
        if (((unsigned int) (str[i] - '0')) > 9)
            return fmpz_set_str_fallback(res, str, base, neg);
    }

    if (slen <= DIGITS_PER_LIMB)
    {
        ulong c = str[0] - '0';

        for (i = 1; i < slen; i++)
            c = c * 10 + (ulong) (str[i] - '0');

        if (neg)
            fmpz_neg_ui(res, c);
        else
            fmpz_set_ui(res, c);
    }
    else
    {
        if (slen < BASECASE_CUTOFF)
            _fmpz_set_str_basecase(res, str, slen);
        else
            fmpz_set_str_bsplit_threaded(res, str, slen);

        if (neg)
            fmpz_neg(res, res);
    }

    return 0;
}
