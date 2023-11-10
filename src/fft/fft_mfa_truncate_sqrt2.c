/*
    Copyright (C) 2009, 2011, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "fft.h"

void fft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v,
    mp_limb_t * s, mp_limb_t * t, mp_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2)
{
   mp_limb_t nw = limbs*FLINT_BITS;
   mp_size_t x, y;
   int negate1 = 0;
   int negate2 = 0;

   if (b1 >= nw)
   {
      negate2 = 1;
      b1 -= nw;
   }
   x  = b1/FLINT_BITS;
   b1 = b1%FLINT_BITS;

   if (b2 >= nw)
   {
      negate1 = 1;
      b2 -= nw;
   }
   y  = b2/FLINT_BITS;
   b2 = b2%FLINT_BITS;

   butterfly_lshB(u, v, s, t, limbs, x, y);
   mpn_mul_2expmod_2expp1(u, u, limbs, b1);
   if (negate2) mpn_neg(u, u, limbs + 1);
   mpn_mul_2expmod_2expp1(v, v, limbs, b2);
   if (negate1) mpn_neg(v, v, limbs + 1);
}

void fft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs)
{
   mp_size_t i;
   mp_size_t limbs = (w*n)/FLINT_BITS;

   if (n == 1)
   {
      mp_size_t tw1 = r*c;
      mp_size_t tw2 = tw1 + rs*c;
      fft_butterfly_twiddle(*t1, *t2, ii[0], ii[is], limbs, tw1*ws, tw2*ws);

      SWAP_PTRS(ii[0],  *t1);
      SWAP_PTRS(ii[is], *t2);

      return;
   }

   for (i = 0; i < n; i++)
   {
      fft_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, limbs, w);

      SWAP_PTRS(ii[i*is],     *t1);
      SWAP_PTRS(ii[(n+i)*is], *t2);
   }

   fft_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, ws, r, c, 2*rs);
   fft_radix2_twiddle(ii+n*is, is, n/2, 2*w, t1, t2, ws, r + rs, c, 2*rs);
}

void fft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)
{
   mp_size_t i;
   mp_size_t limbs = (w*n)/FLINT_BITS;

   if (trunc == 2*n)
      fft_radix2_twiddle(ii, is, n, w, t1, t2, ws, r, c, rs);
   else if (trunc <= n)
   {
      for (i = 0; i < n; i++)
         mpn_add_n(ii[i*is], ii[i*is], ii[(i+n)*is], limbs + 1);

      fft_truncate1_twiddle(ii, is, n/2, 2*w, t1, t2, ws, r, c, 2*rs, trunc);
   } else
   {
      for (i = 0; i < n; i++)
      {
         fft_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, limbs, w);

         SWAP_PTRS(ii[i*is],     *t1);
         SWAP_PTRS(ii[(n+i)*is], *t2);
      }

      fft_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, ws, r, c, 2*rs);
      fft_truncate1_twiddle(ii + n*is, is, n/2, 2*w,
                                     t1, t2, ws, r + rs, c, 2*rs, trunc - n);
   }
}

void fft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n,
                   flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                             mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
{
   mp_size_t i, j, s;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t trunc2 = (trunc - 2*n)/n1;
   mp_size_t limbs = (n*w)/FLINT_BITS;
   flint_bitcnt_t depth = 0;
   flint_bitcnt_t depth2 = 0;

   while ((UWORD(1)<<depth) < n2) depth++;
   while ((UWORD(1)<<depth2) < n1) depth2++;

   /* first half matrix fourier FFT : n2 rows, n1 cols */

   /* FFTs on columns */
   for (i = 0; i < n1; i++)
   {
      /* relevant part of first layer of full sqrt2 FFT */
      if (w & 1)
      {
         for (j = i; j < trunc - 2*n; j+=n1)
         {
            if (j & 1)
               fft_butterfly_sqrt2(*t1, *t2, ii[j], ii[2*n+j], j, limbs, w, *temp);
            else
               fft_butterfly(*t1, *t2, ii[j], ii[2*n+j], j/2, limbs, w);

            SWAP_PTRS(ii[j],     *t1);
            SWAP_PTRS(ii[2*n+j], *t2);
         }

         for ( ; j < 2*n; j+=n1)
         {
             if (i & 1)
                fft_adjust_sqrt2(ii[j + 2*n], ii[j], j, limbs, w, *temp);
             else
                fft_adjust(ii[j + 2*n], ii[j], j/2, limbs, w);
         }
      } else
      {
         for (j = i; j < trunc - 2*n; j+=n1)
         {
            fft_butterfly(*t1, *t2, ii[j], ii[2*n+j], j, limbs, w/2);

            SWAP_PTRS(ii[j],     *t1);
            SWAP_PTRS(ii[2*n+j], *t2);
         }

         for ( ; j < 2*n; j+=n1)
            fft_adjust(ii[j + 2*n], ii[j], j, limbs, w/2);
      }

      /*
         FFT of length n2 on column i, applying z^{r*i} for rows going up in steps
         of 1 starting at row 0, where z => w bits
      */

      fft_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1);
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = n_revbin(j, depth);
         if (j < s) SWAP_PTRS(ii[i+j*n1], ii[i+s*n1]);
      }
   }

   /* FFTs on rows */
   for (i = 0; i < n2; i++)
   {
      fft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = n_revbin(j, depth2);
         if (j < t) SWAP_PTRS(ii[i*n1+j], ii[i*n1+t]);
      }
   }

   /* second half matrix fourier FFT : n2 rows, n1 cols */
   ii += 2*n;

   /* FFTs on columns */
   for (i = 0; i < n1; i++)
   {
      /*
         FFT of length n2 on column i, applying z^{r*i} for rows going up in steps
         of 1 starting at row 0, where z => w bits
      */

      fft_truncate1_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1, trunc2);
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = n_revbin(j, depth);
         if (j < s) SWAP_PTRS(ii[i+j*n1], ii[i+s*n1]);
      }
   }

   /* FFTs on relevant rows */
   for (s = 0; s < trunc2; s++)
   {
      i = n_revbin(s, depth);
      fft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);

      for (j = 0; j < n1; j++)
      {
         mp_size_t t = n_revbin(j, depth2);
         if (j < t) SWAP_PTRS(ii[i*n1+j], ii[i*n1+t]);
      }
   }
}

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
    mp_limb_t ** t1;
    mp_limb_t ** t2;
    mp_limb_t * temp;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
}
fft_outer_arg_t;

void
_fft_outer1_worker(void * arg_ptr)
{
    fft_outer_arg_t arg = *((fft_outer_arg_t *) arg_ptr);
    mp_size_t n1 = arg.n1;
    mp_size_t n2 = arg.n2;
    mp_size_t n = arg.n;
    mp_size_t trunc = arg.trunc;
    mp_size_t limbs = arg.limbs;
    flint_bitcnt_t depth = arg.depth;
    flint_bitcnt_t w = arg.w;
    mp_limb_t ** ii = arg.ii;
    mp_limb_t ** t1 = arg.t1;
    mp_limb_t ** t2 = arg.t2;
    mp_limb_t * temp = arg.temp;
    mp_size_t i, j, end;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
	i = *arg.i;
        end = *arg.i = FLINT_MIN(i + 16, n1);
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (i >= n1)
            return;

        for ( ; i < end; i++)
        {
            /* relevant part of first layer of full sqrt2 FFT */
            if (w & 1)
            {
                for (j = i; j < trunc - 2*n; j+=n1)
                {
                    if (j & 1)
                        fft_butterfly_sqrt2(*t1, *t2, ii[j], ii[2*n+j],
                                                            j, limbs, w, temp);
                    else
                        fft_butterfly(*t1, *t2, ii[j], ii[2*n+j], j/2, limbs, w);

                    SWAP_PTRS(ii[j],     *t1);
                    SWAP_PTRS(ii[2*n+j], *t2);
                }

                for ( ; j < 2*n; j+=n1)
                {
                    if (i & 1)
                        fft_adjust_sqrt2(ii[j + 2*n], ii[j], j, limbs, w, temp);
                    else
                        fft_adjust(ii[j + 2*n], ii[j], j/2, limbs, w);
                }
            } else
            {
                for (j = i; j < trunc - 2*n; j+=n1)
                {
                    fft_butterfly(*t1, *t2, ii[j], ii[2*n+j], j, limbs, w/2);

                    SWAP_PTRS(ii[j],     *t1);
                    SWAP_PTRS(ii[2*n+j], *t2);
                }

                for ( ; j < 2*n; j+=n1)
                    fft_adjust(ii[j + 2*n], ii[j], j, limbs, w/2);
            }

            /*
                FFT of length n2 on column i, applying z^{r*i} for rows going up in steps
                of 1 starting at row 0, where z => w bits
            */

            fft_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1);
            for (j = 0; j < n2; j++)
            {
                mp_size_t s = n_revbin(j, depth);
                if (j < s) SWAP_PTRS(ii[i + j*n1], ii[i + s*n1]);
            }
        }
    }
}

void
_fft_outer2_worker(void * arg_ptr)
{
    fft_outer_arg_t arg = *((fft_outer_arg_t *) arg_ptr);
    mp_size_t n1 = arg.n1;
    mp_size_t n2 = arg.n2;
    mp_size_t trunc2 = arg.trunc;
    flint_bitcnt_t depth = arg.depth;
    flint_bitcnt_t w = arg.w;
    mp_limb_t ** ii = arg.ii;
    mp_limb_t ** t1 = arg.t1;
    mp_limb_t ** t2 = arg.t2;
    mp_size_t i, j, end;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
        i = *arg.i;
        end = *arg.i = FLINT_MIN(i + 16, n1);
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (i >= n1)
            return;

        for ( ; i < end; i++)
        {
            /*
                FFT of length n2 on column i, applying z^{r*i} for rows going up in steps
                of 1 starting at row 0, where z => w bits
            */

            fft_truncate1_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1, trunc2);
            for (j = 0; j < n2; j++)
            {
                mp_size_t s = n_revbin(j, depth);
                if (j < s) SWAP_PTRS(ii[i+j*n1], ii[i+s*n1]);
            }
        }
    }
}

void fft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n,
                   flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                             mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
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
    fft_outer_arg_t * args;

    while ((UWORD(1)<<depth) < n2) depth++;

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    /* first half matrix fourier FFT : n2 rows, n1 cols */

    /* FFTs on columns */

    num_threads = flint_request_threads(&threads,
                             FLINT_MIN(flint_get_num_threads(), (n1 + 15)/16));

    args = (fft_outer_arg_t *)
                       flint_malloc(sizeof(fft_outer_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].i = &shared_i;
       args[i].n1 = n1;
       args[i].n2 = n2;
       args[i].n = n;
       args[i].trunc = trunc;
       args[i].limbs = limbs;
       args[i].depth = depth;
       args[i].w = w;
       args[i].ii = ii;
       args[i].t1 = t1 + i;
       args[i].t2 = t2 + i;
       args[i].temp = temp[i];
#if FLINT_USES_PTHREAD
       args[i].mutex = &mutex;
#endif
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                 _fft_outer1_worker, &args[i]);

    _fft_outer1_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    /* second half matrix fourier FFT : n2 rows, n1 cols */
    ii += 2*n;

    /* FFTs on columns */

    shared_i = 0;


    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].trunc = trunc2;
       args[i].ii = ii;
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                 _fft_outer2_worker, &args[i]);

    _fft_outer2_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

    flint_free(args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif
}
