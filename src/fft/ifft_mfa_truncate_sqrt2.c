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

void ifft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v,
   mp_limb_t * s, mp_limb_t * t, mp_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2)
{
   mp_limb_t nw = limbs*FLINT_BITS;
   mp_size_t x, y;
   int negate1 = 0;
   int negate2 = 0;

   if (b1 >= nw)
   {
      negate1 = 1;
      b1 -= nw;
   }
   x  = b1/FLINT_BITS;
   b1 = b1%FLINT_BITS;

   if (b2 >= nw)
   {
      negate2 = 1;
      b2 -= nw;
   }
   y  = b2/FLINT_BITS;
   b2 = b2%FLINT_BITS;

   if (negate1) mpn_neg(s, s, limbs + 1);
   mpn_div_2expmod_2expp1(s, s, limbs, b1);
   if (negate2) mpn_neg(t, t, limbs + 1);
   mpn_div_2expmod_2expp1(t, t, limbs, b2);
   butterfly_rshB(u, v, s, t, limbs, x, y);
}

void ifft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
        mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                            mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs)
{
   mp_size_t i;
   mp_size_t limbs = (w*n)/FLINT_BITS;

   if (n == 1)
   {
      mp_size_t tw1, tw2;
      tw1 = r*c;
      tw2 = tw1 + rs*c;
      ifft_butterfly_twiddle(*t1, *t2, ii[0], ii[is], limbs, tw1*ws, tw2*ws);

      SWAP_PTRS(ii[0],  *t1);
      SWAP_PTRS(ii[is], *t2);

      return;
   }

   ifft_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, ws, r, c, 2*rs);
   ifft_radix2_twiddle(ii+n*is, is, n/2, 2*w, t1, t2, ws, r + rs, c, 2*rs);

   for (i = 0; i < n; i++)
   {
      ifft_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, limbs, w);

      SWAP_PTRS(ii[i*is], *t1);
      SWAP_PTRS(ii[(n+i)*is], *t2);
   }
}

void ifft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
        mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
           mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)
{
   mp_size_t i;
   mp_size_t limbs = (w*n)/FLINT_BITS;

   if (trunc == 2*n)
      ifft_radix2_twiddle(ii, is, n, w, t1, t2, ws, r, c, rs);
   else if (trunc <= n)
   {
      for (i = trunc; i < n; i++)
      {
         mpn_add_n(ii[i*is], ii[i*is], ii[(i+n)*is], limbs + 1);
         mpn_div_2expmod_2expp1(ii[i*is], ii[i*is], limbs, 1);
      }

      ifft_truncate1_twiddle(ii, is, n/2, 2*w, t1, t2, ws, r, c, 2*rs, trunc);

      for (i = 0; i < trunc; i++)
      {
#if HAVE_ADDSUB_N
         mpn_addsub_n(ii[i*is], ii[i*is], ii[i*is], ii[(n+i)*is], limbs + 1);
#else
         mpn_add_n(ii[i*is], ii[i*is], ii[i*is], limbs + 1);
         mpn_sub_n(ii[i*is], ii[i*is], ii[(n+i)*is], limbs + 1);
#endif
      }
   } else
   {
      ifft_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, ws, r, c, 2*rs);

      for (i = trunc - n; i < n; i++)
      {
          mpn_sub_n(ii[(i+n)*is], ii[i*is], ii[(i+n)*is], limbs + 1);
          fft_adjust(*t1, ii[(i+n)*is], i, limbs, w);
          mpn_add_n(ii[i*is], ii[i*is], ii[(i+n)*is], limbs + 1);
          SWAP_PTRS(ii[(i+n)*is], *t1);
      }

      ifft_truncate1_twiddle(ii + n*is, is, n/2, 2*w, t1, t2, ws, r + rs, c, 2*rs, trunc - n);

      for (i = 0; i < trunc - n; i++)
      {
         ifft_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, limbs, w);

         SWAP_PTRS(ii[i*is],     *t1);
         SWAP_PTRS(ii[(n+i)*is], *t2);
      }
   }
}

void ifft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
   mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
{
   mp_size_t i, j, s;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t trunc2 = (trunc - 2*n)/n1;
   flint_bitcnt_t depth = 0;
   flint_bitcnt_t depth2 = 0;
   flint_bitcnt_t limbs = (w*n)/FLINT_BITS;

   while ((UWORD(1)<<depth) < n2) depth++;
   while ((UWORD(1)<<depth2) < n1) depth2++;

   /* first half mfa IFFT : n2 rows, n1 cols */

   /* row IFFTs */
   for (i = 0; i < n2; i++)
   {
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = n_revbin(j, depth2);
         if (j < s) SWAP_PTRS(ii[i*n1+j], ii[i*n1+s]);
      }

      ifft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);
   }

   /* column IFFTs */
   for (i = 0; i < n1; i++)
   {
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = n_revbin(j, depth);
         if (j < s) SWAP_PTRS(ii[i+j*n1], ii[i+s*n1]);
      }

      /*
         IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps
         of 1 starting at row 0, where z => w bits
      */
      ifft_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1);
   }

   /* second half IFFT : n2 rows, n1 cols */
   ii += 2*n;

   /* row IFFTs */
   for (s = 0; s < trunc2; s++)
   {
      i = n_revbin(s, depth);
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = n_revbin(j, depth2);
         if (j < t) SWAP_PTRS(ii[i*n1 + j], ii[i*n1 + t]);
      }

      ifft_radix2(ii + i*n1, n1/2, w*n2, t1, t2);
   }

   /* column IFFTs with relevant sqrt2 layer butterflies combined */
   for (i = 0; i < n1; i++)
   {
      for (j = 0; j < trunc2; j++)
      {
         mp_size_t s = n_revbin(j, depth);
         if (j < s) SWAP_PTRS(ii[i+j*n1], ii[i+s*n1]);
      }

      for ( ; j < n2; j++)
      {
         mp_size_t u = i + j*n1;
         if (w & 1)
         {
            if (i & 1)
               fft_adjust_sqrt2(ii[i + j*n1], ii[u - 2*n], u, limbs, w, *temp);
            else
               fft_adjust(ii[i + j*n1], ii[u - 2*n], u/2, limbs, w);
         } else
            fft_adjust(ii[i + j*n1], ii[u - 2*n], u, limbs, w/2);
      }

      /*
         IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps
         of 1 starting at row 0, where z => w bits
      */
      ifft_truncate1_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1, trunc2);

      /* relevant components of final sqrt2 layer of IFFT */
      if (w & 1)
      {
         for (j = i; j < trunc - 2*n; j+=n1)
         {
            if (j & 1)
               ifft_butterfly_sqrt2(*t1, *t2, ii[j - 2*n], ii[j], j, limbs, w, *temp);
            else
               ifft_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j/2, limbs, w);

            SWAP_PTRS(ii[j-2*n], *t1);
            SWAP_PTRS(ii[j],     *t2);
         }
      } else
      {
         for (j = i; j < trunc - 2*n; j+=n1)
         {
            ifft_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j, limbs, w/2);

            SWAP_PTRS(ii[j-2*n], *t1);
            SWAP_PTRS(ii[j],     *t2);
         }
      }

      for (j = trunc + i - 2*n; j < 2*n; j+=n1)
           mpn_add_n(ii[j - 2*n], ii[j - 2*n], ii[j - 2*n], limbs + 1);
   }
}

typedef struct
{
    volatile mp_size_t * i;
    mp_size_t n1;
    mp_size_t n2;
    mp_size_t n;
    mp_size_t trunc;
    mp_size_t trunc2;
    mp_size_t limbs;
    flint_bitcnt_t depth;
    flint_bitcnt_t depth2;
    flint_bitcnt_t w;
    mp_limb_t ** ii;
    mp_limb_t ** t1;
    mp_limb_t ** t2;
    mp_limb_t * temp;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
}
ifft_outer_arg_t;

void
_ifft_outer1_worker(void * arg_ptr)
{
    ifft_outer_arg_t arg = *((ifft_outer_arg_t *) arg_ptr);
    mp_size_t n1 = arg.n1;
    mp_size_t n2 = arg.n2;
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
            for (j = 0; j < n2; j++)
            {
                mp_size_t s = n_revbin(j, depth);
                if (j < s) SWAP_PTRS(ii[i + j*n1], ii[i + s*n1]);
            }

            /*
                IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps
                of 1 starting at row 0, where z => w bits
            */
            ifft_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1);
        }
    }
}

void
_ifft_outer2_worker(void * arg_ptr)
{
    ifft_outer_arg_t arg = *((ifft_outer_arg_t *) arg_ptr);
    mp_size_t n1 = arg.n1;
    mp_size_t n2 = arg.n2;
    mp_size_t n = arg.n;
    mp_size_t trunc = arg.trunc;
    mp_size_t trunc2 = arg.trunc2;
    mp_size_t limbs = arg.limbs;
    flint_bitcnt_t depth = arg.depth;
    flint_bitcnt_t depth2 = arg.depth2;
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
            for (j = 0; j < trunc2; j++)
            {
                mp_size_t s = n_revbin(j, depth);
                if (j < s) SWAP_PTRS(ii[i + j*n1], ii[i + s*n1]);
            }

            for ( ; j < n2; j++)
            {
                mp_size_t u = i + j*n1;
                if (w & 1)
                {
                    if (i & 1)
                        fft_adjust_sqrt2(ii[i + j*n1], ii[u - 2*n], u, limbs, w, temp);
                    else
                        fft_adjust(ii[i + j*n1], ii[u - 2*n], u/2, limbs, w);
                } else
                    fft_adjust(ii[i + j*n1], ii[u - 2*n], u, limbs, w/2);
            }

            /*
                IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps
                of 1 starting at row 0, where z => w bits
            */
            ifft_truncate1_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, w, 0, i, 1, trunc2);

            /* relevant components of final sqrt2 layer of IFFT */
            if (w & 1)
            {
                for (j = i; j < trunc - 2*n; j+=n1)
                {
                    if (j & 1)
                        ifft_butterfly_sqrt2(*t1, *t2, ii[j - 2*n], ii[j], j, limbs, w, temp);
                    else
                        ifft_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j/2, limbs, w);

                    SWAP_PTRS(ii[j - 2*n], *t1);
                    SWAP_PTRS(ii[j],       *t2);
                }
            } else
            {
                for (j = i; j < trunc - 2*n; j+=n1)
                {
                    ifft_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j, limbs, w/2);

                    SWAP_PTRS(ii[j - 2*n], *t1);
                    SWAP_PTRS(ii[j],       *t2);
                }
            }

            for (j = trunc + i - 2*n; j < 2*n; j+=n1)
                mpn_add_n(ii[j - 2*n], ii[j - 2*n], ii[j - 2*n], limbs + 1);

            for (j = 0; j < trunc2; j++)
            {
                mp_size_t t = j*n1 + i;
                mpn_div_2expmod_2expp1(ii[t], ii[t], limbs, depth + depth2 + 1);
                mpn_normmod_2expp1(ii[t], limbs);
            }

            for (j = 0; j < n2; j++)
            {
                mp_size_t t = j*n1 + i - 2*n;
                mpn_div_2expmod_2expp1(ii[t], ii[t], limbs, depth + depth2 + 1);
                mpn_normmod_2expp1(ii[t], limbs);
            }
        }
    }
}

void ifft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
   mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
{
    mp_size_t i, shared_i = 0;
    mp_size_t n2 = (2*n)/n1;
    mp_size_t trunc2 = (trunc - 2*n)/n1;
    flint_bitcnt_t depth = 0;
    flint_bitcnt_t depth2 = 0;
    flint_bitcnt_t limbs = (w*n)/FLINT_BITS;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong num_threads;
    thread_pool_handle * threads;
    ifft_outer_arg_t * args;

    while ((UWORD(1)<<depth) < n2) depth++;
    while ((UWORD(1)<<depth2) < n1) depth2++;

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

   /* first half mfa IFFT : n2 rows, n1 cols */

   /* column IFFTs */

    num_threads = flint_request_threads(&threads,
                             FLINT_MIN(flint_get_num_threads(), (n1 + 15)/16));

    args = (ifft_outer_arg_t *)
                      flint_malloc(sizeof(ifft_outer_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].i = &shared_i;
       args[i].n1 = n1;
       args[i].n2 = n2;
       args[i].n = n;
       args[i].trunc = trunc;
       args[i].trunc2 = trunc2;
       args[i].limbs = limbs;
       args[i].depth = depth;
       args[i].depth2 = depth2;
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
                                                _ifft_outer1_worker, &args[i]);

    _ifft_outer1_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    /* second half IFFT : n2 rows, n1 cols */
    ii += 2*n;

    /* column IFFTs with relevant sqrt2 layer butterflies combined */

    shared_i = 0;

    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].ii = ii;
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                _ifft_outer2_worker, &args[i]);

    _ifft_outer2_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

    flint_free(args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif
}
