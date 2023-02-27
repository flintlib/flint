/*
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <pthread.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "thread_support.h"

typedef struct
{
    fmpz * vec;
    mp_ptr * residues;
    slong n0;
    slong n1;
    mp_srcptr primes;
    slong num_primes;
    int crt;  /* reduce if 0, lift if 1 */
}
mod_ui_arg_t;

void
_fmpz_vec_multi_mod_ui_worker(void * arg_ptr)
{
    mod_ui_arg_t arg = *((mod_ui_arg_t *) arg_ptr);
    mp_ptr tmp;
    slong i, j;

    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;

    tmp = flint_malloc(sizeof(mp_limb_t) * arg.num_primes);
    fmpz_comb_init(comb, arg.primes, arg.num_primes);
    fmpz_comb_temp_init(comb_temp, comb);

    for (i = arg.n0; i < arg.n1; i++)
    {
        if (arg.crt)
        {
            for (j = 0; j < arg.num_primes; j++)
                tmp[j] = arg.residues[j][i];
            fmpz_multi_CRT_ui(arg.vec + i, tmp, comb, comb_temp, 1);
        }
        else
        {
            fmpz_multi_mod_ui(tmp, arg.vec + i, comb, comb_temp);
            for (j = 0; j < arg.num_primes; j++)
                arg.residues[j][i] = tmp[j];
        }
    }

    flint_free(tmp);
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(comb_temp);
}

void
_fmpz_vec_multi_mod_ui_threaded(mp_ptr * residues, fmpz * vec, slong len,
    mp_srcptr primes, slong num_primes, int crt)
{
    mod_ui_arg_t * args;
    slong i, num_threads;
    thread_pool_handle * threads;

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    args = (mod_ui_arg_t *)
                          flint_malloc(sizeof(mod_ui_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].vec = vec;
        args[i].residues = residues;
        args[i].n0 = (len * i) / (num_threads + 1);
        args[i].n1 = (len * (i + 1)) / (num_threads + 1);
        args[i].primes = (mp_ptr) primes;
        args[i].num_primes = num_primes;
        args[i].crt = crt;
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                               _fmpz_vec_multi_mod_ui_worker, &args[i]);

    _fmpz_vec_multi_mod_ui_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

    flint_free(args);
}

typedef struct
{
    mp_ptr * residues;
    slong len;
    mp_srcptr primes;
    slong num_primes;
    slong p0;
    slong p1;
    fmpz * c;
}
taylor_shift_arg_t;

void
_fmpz_poly_multi_taylor_shift_worker(void * arg_ptr)
{
    taylor_shift_arg_t arg = *((taylor_shift_arg_t *) arg_ptr);
    slong i;

    for (i = arg.p0; i < arg.p1; i++)
    {
        nmod_t mod;
        mp_limb_t p, cm;

        p = arg.primes[i];
        nmod_init(&mod, p);
        cm = fmpz_fdiv_ui(arg.c, p);
        _nmod_poly_taylor_shift(arg.residues[i], cm, arg.len, mod);
    }
}

void
_fmpz_poly_multi_taylor_shift_threaded(mp_ptr * residues, slong len,
        const fmpz_t c, mp_srcptr primes, slong num_primes)
{
    taylor_shift_arg_t * args;
    slong i, num_threads;
    thread_pool_handle * threads;

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    args = (taylor_shift_arg_t *)
                    flint_malloc(sizeof(taylor_shift_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].residues = residues;
        args[i].len = len;
        args[i].p0 = (num_primes * i) / (num_threads + 1);
        args[i].p1 = (num_primes * (i + 1)) / (num_threads + 1);
        args[i].primes = (mp_ptr) primes;
        args[i].num_primes = num_primes;
        args[i].c = (fmpz *) c;
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                               _fmpz_poly_multi_taylor_shift_worker, &args[i]);

    _fmpz_poly_multi_taylor_shift_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

    flint_free(args);
}

void
_fmpz_poly_taylor_shift_multi_mod(fmpz * poly, const fmpz_t c, slong len)
{
    slong xbits, ybits, num_primes, i;
    mp_ptr primes;
    mp_ptr * residues;

    if (len <= 1 || fmpz_is_zero(c))
        return;

    xbits = _fmpz_vec_max_bits(poly, len);

    if (xbits == 0)
        return;

    /* If poly has degree D and coefficients at most |C|, the
       output has coefficient at most D * |C| * 2^D * c^D */
    xbits = FLINT_ABS(xbits) + 1;
    ybits = xbits + len + FLINT_BIT_COUNT(len);

    if (!fmpz_is_pm1(c))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_pow_ui(t, c, len);
        ybits += fmpz_bits(t);
        fmpz_clear(t);
    }

    /* Use primes greater than 2^(FLINT_BITS-1) */
    num_primes = (ybits + (FLINT_BITS - 1) - 1) / (FLINT_BITS - 1);
    primes = flint_malloc(sizeof(mp_limb_t) * num_primes);
    primes[0] = n_nextprime(UWORD(1) << (FLINT_BITS - 1), 1);
    for (i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i-1], 1);

    /* Space for poly reduced modulo the primes */
    residues = flint_malloc(sizeof(mp_ptr) * num_primes);
    for (i = 0; i < num_primes; i++)
        residues[i] = flint_malloc(sizeof(mp_limb_t) * len);

    _fmpz_vec_multi_mod_ui_threaded(residues, poly, len, primes,
                                                                num_primes, 0);
    _fmpz_poly_multi_taylor_shift_threaded(residues, len, c,
                                                           primes, num_primes);
    _fmpz_vec_multi_mod_ui_threaded(residues, poly, len, primes,
		                                                num_primes, 1);

    for (i = 0; i < num_primes; i++)
        flint_free(residues[i]);
    flint_free(residues);
    flint_free(primes);
}
