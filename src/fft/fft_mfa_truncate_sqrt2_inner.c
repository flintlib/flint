/*
    Copyright (C) 2009, 2011, 2020 William Hart
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"
#include "thread_pool.h"

struct sd_fft_mpn_mulmod_2expp1_ctx_struct;
struct sd_fft_mpn_mulmod_2expp1_scratch_struct;
#include "thread_support.h"
#include "ulong_extras.h"
#include "fft.h"

typedef struct
{
    volatile mp_size_t * i;
    mp_size_t n1;
    mp_size_t n2;
    mp_size_t n;
    mp_size_t trunc;
    mp_size_t limbs;
    flint_bitcnt_t depth;
    flint_bitcnt_t w;
    mp_limb_t ** ii;
    mp_limb_t ** jj;
    mp_limb_t ** t1;
    mp_limb_t ** t2;
    mp_limb_t * tt;
    struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx;
    struct sd_fft_mpn_mulmod_2expp1_scratch_struct * nscr;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
}
fft_inner_arg_t;

static void
_fft_inner1_worker(void * arg_ptr)
{
    fft_inner_arg_t arg = *((fft_inner_arg_t *) arg_ptr);
    mp_size_t n1 = arg.n1;
    mp_size_t n2 = arg.n2;
    mp_size_t n = arg.n;
    mp_size_t trunc = arg.trunc;
    mp_size_t limbs = arg.limbs;
    flint_bitcnt_t depth = arg.depth;
    flint_bitcnt_t w = arg.w;
    mp_limb_t ** ii = arg.ii;
    mp_limb_t ** jj = arg.jj;
    mp_limb_t ** t1 = arg.t1;
    mp_limb_t ** t2 = arg.t2;
    mp_limb_t * tt = arg.tt;
    mp_size_t i, j, s, end;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
        s = *arg.i;
        end = *arg.i = FLINT_MIN(s + 16, trunc);
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (s >= trunc)
            return;

        for ( ; s < end; s++)
        {
            i = n_revbin(s, depth);
            fft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);
            if (ii != jj) fft_radix2(jj + i*n1, n1/2, w*n2, t1, t2);

            for (j = 0; j < n1; j++)
            {
                mp_size_t t = i*n1 + j;
                mpn_normmod_2expp1(ii[t], limbs);
                if (ii != jj) mpn_normmod_2expp1(jj[t], limbs);
                #if FLINT_HAVE_FFT_SMALL
                if (arg.sdctx != NULL)
                    sd_fft_mpn_mulmod_2expp1(arg.sdctx, ii[t], ii[t], jj[t],
                                arg.nscr);
                else
#endif
                    fft_mulmod_2expp1(ii[t], ii[t], jj[t], n, w, tt);
            }

            ifft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);
        }
    }
}

static void
_fft_inner2_worker(void * arg_ptr)
{
    fft_inner_arg_t arg = *((fft_inner_arg_t *) arg_ptr);
    mp_size_t n1 = arg.n1;
    mp_size_t n2 = arg.n2;
    mp_size_t n = arg.n;
    mp_size_t limbs = arg.limbs;
    flint_bitcnt_t w = arg.w;
    mp_limb_t ** ii = arg.ii;
    mp_limb_t ** jj = arg.jj;
    mp_limb_t ** t1 = arg.t1;
    mp_limb_t ** t2 = arg.t2;
    mp_limb_t * tt = arg.tt;
    mp_size_t i, j, end;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
        i = *arg.i;
        end = *arg.i = FLINT_MIN(i + 16, n2);
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (i >= n2)
            return;

        for ( ; i < end; i++)
        {
            fft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);
            if (ii != jj) fft_radix2(jj + i*n1, n1/2, w*n2, t1, t2);

            for (j = 0; j < n1; j++)
            {
                mp_size_t t = i*n1 + j;
                mpn_normmod_2expp1(ii[t], limbs);
                if (ii != jj) mpn_normmod_2expp1(jj[t], limbs);
                #if FLINT_HAVE_FFT_SMALL
                if (arg.sdctx != NULL)
                    sd_fft_mpn_mulmod_2expp1(arg.sdctx, ii[t], ii[t], jj[t],
                                arg.nscr);
                else
#endif
                    fft_mulmod_2expp1(ii[t], ii[t], jj[t], n, w, tt);
            }

            ifft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);
        }
    }
}

static void
_fft_mfa_truncate_sqrt2_inner_impl(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx,
                   mp_limb_t ** ii, mp_limb_t ** jj, mp_size_t n,
                   flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                  mp_limb_t ** FLINT_UNUSED(temp), mp_size_t n1, mp_size_t trunc, mp_limb_t ** tt)
{
    mp_size_t i, shared_i = 0;
    mp_size_t n2 = (2*n)/n1;
    mp_size_t trunc2 = (trunc - 2*n)/n1;
    mp_size_t limbs = (n*w)/FLINT_BITS;
    flint_bitcnt_t depth = 0;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong num_threads;
    thread_pool_handle * threads;
    fft_inner_arg_t * args;

    while ((UWORD(1)<<depth) < (ulong) n2) depth++;

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    ii += 2*n;
    jj += 2*n;

   /* convolutions on relevant rows */

    /* Phase 1 has 2n/n1 columns and phase 2 has trunc2 rows of
       work-stolen units; each unit is heavy at large limb counts
       (n1 column transforms plus n1 pointwise multiplications), so
       for the negacyclic engine one unit per thread is worthwhile.
       The classic cap is kept unchanged. */
    num_threads = flint_request_threads(&threads,
                     sdctx != NULL
                       ? FLINT_MIN(flint_get_num_threads(),
                                   FLINT_MAX((2*n)/n1, trunc2))
                       : FLINT_MIN(flint_get_num_threads(),
                                   (trunc2 + 15)/16));

    args = (fft_inner_arg_t *)
                       flint_malloc(sizeof(fft_inner_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].i = &shared_i;
       args[i].n1 = n1;
       args[i].n2 = n2;
       args[i].n = n;
       args[i].trunc = trunc2;
       args[i].limbs = limbs;
       args[i].depth = depth;
       args[i].w = w;
       args[i].ii = ii;
       args[i].jj = jj;
       args[i].t1 = t1 + i;
       args[i].t2 = t2 + i;
       args[i].tt = tt[i];
       args[i].sdctx = sdctx;
       args[i].nscr = NULL;
#if FLINT_USES_PTHREAD
       args[i].mutex = &mutex;
#endif
    }

#if FLINT_HAVE_FFT_SMALL
    if (sdctx != NULL)
    {
        for (i = 0; i < num_threads + 1; i++)
        {
            args[i].nscr = flint_malloc(
                             sizeof(sd_fft_mpn_mulmod_2expp1_scratch_struct));
            sd_fft_mpn_mulmod_2expp1_scratch_init(args[i].nscr, sdctx);
        }
    }
#endif

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                 _fft_inner1_worker, &args[i]);

    _fft_inner1_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    ii -= 2*n;
    jj -= 2*n;

    /* convolutions on rows */

    shared_i = 0;

    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].ii = ii;
       args[i].jj = jj;
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                 _fft_inner2_worker, &args[i]);

    _fft_inner2_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

#if FLINT_HAVE_FFT_SMALL
    if (sdctx != NULL)
    {
        for (i = 0; i < num_threads + 1; i++)
        {
            sd_fft_mpn_mulmod_2expp1_scratch_clear(args[i].nscr);
            flint_free(args[i].nscr);
        }
    }
#endif
    flint_free(args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif
}

void fft_mfa_truncate_sqrt2_inner(mp_limb_t ** ii, mp_limb_t ** jj, mp_size_t n,
                   flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1,
                                mp_size_t n1, mp_size_t trunc, mp_limb_t ** tt)
{
    _fft_mfa_truncate_sqrt2_inner_impl(NULL, ii, jj, n, w, t1, t2, s1,
                                       n1, trunc, tt);
}

#if FLINT_HAVE_FFT_SMALL
/* as fft_mfa_truncate_sqrt2_inner with the pointwise multiplications
   performed by the negacyclic engine; each worker gets its own
   scratch, the context being read-only during multiplication */
void fft_mfa_truncate_sqrt2_inner_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx,
                   mp_limb_t ** ii, mp_limb_t ** jj, mp_size_t n,
                   flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1,
                                mp_size_t n1, mp_size_t trunc, mp_limb_t ** tt)
{
    _fft_mfa_truncate_sqrt2_inner_impl(sdctx, ii, jj, n, w, t1, t2, s1,
                                       n1, trunc, tt);
}
#endif
