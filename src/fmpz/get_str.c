/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "string.h"

/*
    Notes:

    * Right now the purpose of the recursive fmpz code is just to have
      something that takes advantage of the small FFT and is
      multithreaded for huge (million-digit) operands.

    * This code should be rewritten mpn style with an efficient
      mpn basecase. GMP handles the basecase just fine, of course, but
      it's wasteful to let GMP recompute the powers on each basecase
      call. (There are also some tricks that GMP doesn't use.)

    * Incidentally, because fmpzs are normalized, this code
      performs better when the input has long strings of 0s.
      (Fun exercise: make it equally fast for long strings of 9s too.)

    * Other potential improvements: cache transforms, trim the low
      zeros in 10^N when dividing.

    * It is mildly annoying that the output is big endian,
      as we have to shift the string when fmpz_sizeinbase is too big.
      This shift can avoided if we always compute the high part
      first, but that gets tricker when we want multithreading.

    * This is currently base 10 only. Trivial to generalize,
      but who cares about other bases?
*/

/* Falling back to GMP. */
#define STR_BASECASE_CUTOFF_BITS 30000

/* We compute inverses of 10^N to speed up divisions, but it is
   pointless to do so near the root of the tree where each
   power is used only a couple of times. */
#define PREINV_DEPTH 3


typedef struct
{
    char * s;
    slong num_digits;
    fmpz * f;
    const slong * exps;
    slong cur_depth;
    slong depth;
    const fmpz * pows;
    const fmpz_preinvn_struct * preinv;
}
worker_args_struct;

static void
_fmpz_get_str_recursive(char * s, slong num_digits, const fmpz_t f,
    const slong * exps, slong cur_depth, slong depth,
    const fmpz * pows, const fmpz_preinvn_struct * preinv);

static void
worker(void * arg)
{
    worker_args_struct * X = (worker_args_struct * ) arg;
    _fmpz_get_str_recursive(X->s, X->num_digits, X->f, X->exps, X->cur_depth, X->depth, X->pows, X->preinv);
    fmpz_clear(X->f);
}

static void
_fmpz_get_str_recursive(char * s, slong num_digits, const fmpz_t f, const slong * exps, slong cur_depth, slong depth, const fmpz * pows, const fmpz_preinvn_struct * preinv)
{
    if (cur_depth >= depth || fmpz_bits(f) < STR_BASECASE_CUTOFF_BITS)
    {
        char * tmp;
        slong n;
        tmp = fmpz_get_str(NULL, 10, f);
        n = strlen(tmp);
        memcpy(s + num_digits - n, tmp, n);
        flint_free(tmp);
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

        if (cur_depth >= PREINV_DEPTH)
            fmpz_fdiv_qr_preinvn(q, r, f, pows + cur_depth, preinv + cur_depth);
        else
            fmpz_fdiv_qr(q, r, f, pows + cur_depth);

        low_digits->s = s + num_digits - num_right;
        low_digits->num_digits = num_right;
        low_digits->f = r;
        low_digits->exps = exps;
        low_digits->cur_depth = cur_depth + 1;
        low_digits->depth = depth;
        low_digits->pows = pows;
        low_digits->preinv = preinv;

        high_digits->s = s;
        high_digits->num_digits = num_digits - num_right;
        high_digits->f = q;
        high_digits->exps = exps;
        high_digits->cur_depth = cur_depth + 1;
        high_digits->depth = depth;
        high_digits->pows = pows;
        high_digits->preinv = preinv;

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
    }
}

char *
fmpz_get_str_bsplit_threaded(char * s, const fmpz_t f)
{
    slong n, k, depth, leading_zeros;
    slong exps[FLINT_BITS];
    fmpz * pows;
    fmpz_preinvn_struct * preinv;

    if (!COEFF_IS_MPZ(*f))
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (s == NULL)
        s = flint_malloc(mpz_sizeinbase(COEFF_TO_PTR(*f), 10) + 2);

    if (fmpz_sgn(f) < 0)
    {
        /* make a shallow mpz copy */
        fmpz_t g;
        mpz_t u;
        *u = *COEFF_TO_PTR(*f);
        *g = PTR_TO_COEFF(u);
        u->_mp_size = -u->_mp_size;
        s[0] = '-';
        fmpz_get_str_bsplit_threaded(s + 1, g);
        return s;
    }

    n = fmpz_sizeinbase(f, 10);

    exps[0] = (n + 1) / 2;
    depth = 1;

    while (exps[depth - 1] > 2 * 0.301 * STR_BASECASE_CUTOFF_BITS)
    {
        exps[depth] = (exps[depth - 1] + 1) / 2;
        depth++;
    }

    pows = _fmpz_vec_init(depth);
    preinv = flint_malloc(sizeof(fmpz_preinvn_struct) * depth);

    fmpz_ui_pow_ui(pows + depth - 1, 5, exps[depth - 1]);
    for (k = depth - 2; k >= 0; k--)
    {
        fmpz_mul(pows + k, pows + k + 1, pows + k + 1);
        if (exps[k] != 2 * exps[k + 1])
            fmpz_divexact_ui(pows + k, pows + k, 5);
    }

    for (k = 0; k < depth; k++)
        fmpz_mul_2exp(pows + k, pows + k, exps[k]);

    for (k = depth - 1; k >= PREINV_DEPTH; k--)
        fmpz_preinvn_init(preinv + k, pows + k);

    memset(s, '0', n);

    _fmpz_get_str_recursive(s, n, f, exps, 0, depth, pows, preinv);
    leading_zeros = 0;
    while (s[leading_zeros] == '0')
        leading_zeros++;

    if (leading_zeros != 0)
    {
        n -= leading_zeros;
        for (k = 0; k < n; k++)
            s[k] = s[k + leading_zeros];
    }

    s[n] = '\0';

    for (k = depth - 1; k >= PREINV_DEPTH; k--)
        fmpz_preinvn_clear(preinv + k);

    _fmpz_vec_clear(pows, depth);
    flint_free(preinv);

    return s;
}

char * fmpz_get_str(char * str, int b, const fmpz_t f)
{
    FLINT_ASSERT(b >= 2 && b <= 62);

    if (!COEFF_IS_MPZ(*f))
    {
        fmpz c;
        mp_limb_t d;
        c = *f;

        d = FLINT_ABS(c);

        /* Need a special case for zero, which may as well handle
           single digits. */
        if (d < FLINT_MIN(b, 10))
        {
            if (str == NULL)
                str = flint_malloc(3);
            str[0] = '-';
            str[0 + (c < 0)] = d + '0';
            str[1 + (c < 0)] = '\0';
        }
        else if (b == 10)
        {
            unsigned char tmp[FLINT_BITS + 3];
            slong i, len;
            /* The compiler might generate faster code for 32-bit divisions */
            unsigned int dl;

            len = 0;
#if FLINT_BITS == 64
            while (d >= (UWORD(1) << 32))
            {
                tmp[len] = d % 10;
                d /= 10;
                len++;
            }
#endif
            dl = d;
            while (dl != 0)
            {
                tmp[len] = dl % 10;
                dl /= 10;
                len++;
            }

            if (str == NULL)
                str = flint_malloc(len + 2);

            str[0] = '-';
            for (i = 0; i < len; i++)
                str[i + (c < 0)] = tmp[len - 1 - i] + '0';
            str[len + (c < 0)] = '\0';

            return str;
        }
        else
        {
            mpz_t z;

            z->_mp_d = &d;
            z->_mp_alloc = 1;
            z->_mp_size = (c > 1) ? 1 : -1;

            if (str == NULL)
                str = flint_malloc(mpz_sizeinbase(z, b) + 2);
            str = mpz_get_str(str, b, z);
        }
    }
    else
    {
        if (str == NULL)
            str = flint_malloc(mpz_sizeinbase(COEFF_TO_PTR(*f), b) + 2);

#ifdef FLINT_HAVE_FFT_SMALL
        if (b == 10 && mpz_size(COEFF_TO_PTR(*f)) > 15000)
        {
            fmpz_get_str_bsplit_threaded(str, f);
        }
        else
#endif
        {
            str = mpz_get_str(str, b, COEFF_TO_PTR(*f));
        }
    }

    return str;
}
