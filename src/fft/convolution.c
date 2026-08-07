/*
    Copyright (C) 2008-2011, 2020 William Hart
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fft.h"
#include "fft_small.h"
#include "thread_pool.h"
#include "thread_support.h"

void fft_convolution(mp_limb_t ** ii, mp_limb_t ** jj, slong depth,
                              slong limbs, slong trunc, mp_limb_t ** t1,
                          mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)
{
   slong n = (WORD(1)<<depth), j;
   slong w = (limbs*FLINT_BITS)/n;
   slong sqrt = (WORD(1)<<(depth/2));

   if (depth <= 6)
   {
      trunc = 2*((trunc + 1)/2);

      fft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      if (ii != jj)
         fft_truncate_sqrt2(jj, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         if (ii != jj) mpn_normmod_2expp1(jj[j], limbs);

         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
      }

      ifft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
         mpn_normmod_2expp1(ii[j], limbs);
      }
   } else
   {
      trunc = 2*sqrt*((trunc + 2*sqrt - 1)/(2*sqrt));

      fft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);

      if (ii != jj)
         fft_mfa_truncate_sqrt2_outer(jj, n, w, t1, t2, s1, sqrt, trunc);

      fft_mfa_truncate_sqrt2_inner(ii, jj, n, w, t1, t2, s1, sqrt, trunc, tt);

      ifft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);
   }
}

/* less optimal version of the above, for test code only */
void fft_convolution_basic(mp_limb_t ** ii, mp_limb_t ** jj, slong depth,
                              slong limbs, slong trunc, mp_limb_t ** t1,
                          mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)
{
   slong n = (WORD(1)<<depth), j, s, t, u, trunc2;
   slong w = (limbs*FLINT_BITS)/n;
   slong sqrt = (WORD(1)<<(depth/2));

   if (depth <= 6)
   {
      trunc = 2*((trunc + 1)/2);

      fft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      if (ii != jj)
         fft_truncate_sqrt2(jj, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         if (ii != jj) mpn_normmod_2expp1(jj[j], limbs);

         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
      }

      ifft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
         mpn_normmod_2expp1(ii[j], limbs);
      }
   } else
   {
      trunc = 2*sqrt*((trunc + 2*sqrt - 1)/(2*sqrt));

      fft_mfa_truncate_sqrt2(ii, n, w, t1, t2, s1, sqrt, trunc);

      if (ii != jj)
         fft_mfa_truncate_sqrt2(jj, n, w, t1, t2, s1, sqrt, trunc);

      for (j = 0; j < 2*n; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         if (ii != jj) mpn_normmod_2expp1(jj[j], limbs);

         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
      }

      trunc2 = (trunc - 2*n)/sqrt;

      for (s = 0; s < trunc2; s++)
      {
         t = n_revbin(s, depth - depth/2 + 1);

         for (u = 0; u < sqrt; u++)
         {
            j = 2*n + t*sqrt + u;
            mpn_normmod_2expp1(ii[j], limbs);
            if (ii != jj) mpn_normmod_2expp1(jj[j], limbs);

            fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
         }
      }

      ifft_mfa_truncate_sqrt2(ii, n, w, t1, t2, s1, sqrt, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
         mpn_normmod_2expp1(ii[j], limbs);
      }
   }
}

#if FLINT_HAVE_FFT_SMALL
/* Threaded pointwise stage of the negacyclic convolutions. The
   coefficient work is embarrassingly parallel: the context is
   read-only during multiplication and each worker gets its own
   scratch. Two range shapes are needed: a contiguous [j0, j1), and
   the truncated-MFA row set j = base + n_revbin(s, rvbits)*srt + u
   for s in [s0, s1), u in [0, srt). */
typedef struct
{
    sd_fft_mpn_mulmod_2expp1_ctx_struct * C;
    mp_limb_t ** ii;
    mp_limb_t ** jj;
    slong j0;
    slong j1;
    slong base;
    slong srt;
    slong rvbits;
    slong limbs;
    sd_fft_mpn_mulmod_2expp1_scratch_struct * scr;
    int rows;
}
_sd_fft_pw_arg_struct;

static void
_sd_fft_pw_one(_sd_fft_pw_arg_struct * X, slong j)
{
    mpn_normmod_2expp1(X->ii[j], X->limbs);
    if (X->ii != X->jj)
        mpn_normmod_2expp1(X->jj[j], X->limbs);
    sd_fft_mpn_mulmod_2expp1(X->C, X->ii[j], X->ii[j], X->jj[j], X->scr);
}

static void
_sd_fft_pw_worker(void * varg)
{
    _sd_fft_pw_arg_struct * X = (_sd_fft_pw_arg_struct *) varg;
    slong s, u, j;

    if (!X->rows)
    {
        for (j = X->j0; j < X->j1; j++)
            _sd_fft_pw_one(X, j);
    }
    else
    {
        for (s = X->j0; s < X->j1; s++)
        {
            slong t = (slong) n_revbin((ulong) s, (ulong) X->rvbits);
            for (u = 0; u < X->srt; u++)
                _sd_fft_pw_one(X, X->base + t*X->srt + u);
        }
    }
}

static void
_sd_fft_pw_run(_sd_fft_pw_arg_struct * proto, slong count)
{
    thread_pool_handle * handles;
    slong nhandles, nthreads, i;
    _sd_fft_pw_arg_struct * args;
    sd_fft_mpn_mulmod_2expp1_scratch_struct * scr;
    slong j0 = proto->j0;

    if (count <= 0)
        return;

    nhandles = flint_request_threads(&handles, FLINT_MIN(count,
                                     flint_get_num_threads()));
    nthreads = nhandles + 1;

    if (nhandles == 0)
    {
        sd_fft_mpn_mulmod_2expp1_scratch_t S;
        sd_fft_mpn_mulmod_2expp1_scratch_init(S, proto->C);
        proto->scr = S;
        _sd_fft_pw_worker(proto);
        sd_fft_mpn_mulmod_2expp1_scratch_clear(S);
        flint_give_back_threads(handles, nhandles);
        return;
    }

    args = flint_malloc(nthreads * sizeof(_sd_fft_pw_arg_struct));
    scr = flint_malloc(nthreads
              * sizeof(sd_fft_mpn_mulmod_2expp1_scratch_struct));
    for (i = 0; i < nthreads; i++)
        sd_fft_mpn_mulmod_2expp1_scratch_init(scr + i, proto->C);

    for (i = 0; i < nthreads; i++)
    {
        args[i] = *proto;
        args[i].j0 = j0 + (i * count) / nthreads;
        args[i].j1 = j0 + ((i + 1) * count) / nthreads;
        args[i].scr = scr + i;
    }

    for (i = nhandles; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _sd_fft_pw_worker, args + i);
    _sd_fft_pw_worker(args + 0);
    for (i = nhandles; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    for (i = 0; i < nthreads; i++)
        sd_fft_mpn_mulmod_2expp1_scratch_clear(scr + i);
    flint_free(scr);
    flint_free(args);
    flint_give_back_threads(handles, nhandles);
}

void
_fft_pointwise_sd_fft(sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii,
                      mp_limb_t ** jj, slong j0, slong j1, slong limbs)
{
    _sd_fft_pw_arg_struct X = { sdctx, ii, jj, j0, j1,
                                0, 0, 0, limbs, NULL, 0 };
    _sd_fft_pw_run(&X, j1 - j0);
}

void
_fft_pointwise_rows_sd_fft(sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx,
                      mp_limb_t ** ii, mp_limb_t ** jj, slong base,
                      slong srt, slong nrows, slong rvbits, slong limbs)
{
    _sd_fft_pw_arg_struct X = { sdctx, ii, jj, 0, nrows,
                                base, srt, rvbits, limbs, NULL, 1 };
    _sd_fft_pw_run(&X, nrows);
}
#else
/* without the fft_small engine no context can exist; these are
   never reached and exist only to satisfy references */
void
_fft_pointwise_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx,
                      mp_limb_t ** ii, mp_limb_t ** jj, slong j0,
                      slong j1, slong limbs)
{
    FLINT_ASSERT(sdctx == NULL);
    (void) sdctx; (void) ii; (void) jj;
    (void) j0; (void) j1; (void) limbs;
}

void
_fft_pointwise_rows_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx,
                      mp_limb_t ** ii, mp_limb_t ** jj, slong base,
                      slong srt, slong nrows, slong rvbits, slong limbs)
{
    FLINT_ASSERT(sdctx == NULL);
    (void) sdctx; (void) ii; (void) jj; (void) base;
    (void) srt; (void) nrows; (void) rvbits; (void) limbs;
}
#endif

/* fft_convolution_basic with an optional negacyclic pointwise
   engine; sdctx == NULL reproduces the classic behavior */
void fft_convolution_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx,
                          mp_limb_t ** ii, mp_limb_t ** jj, slong depth,
                              slong limbs, slong trunc, mp_limb_t ** t1,
                          mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)
{
   slong n = (WORD(1)<<depth), j;
   slong w = (limbs*FLINT_BITS)/n;
   slong sqrt = (WORD(1)<<(depth/2));

   if (depth <= 6)
   {
      trunc = 2*((trunc + 1)/2);

      fft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      if (ii != jj)
         fft_truncate_sqrt2(jj, n, w, t1, t2, s1, trunc);

#if FLINT_HAVE_FFT_SMALL
      if (sdctx != NULL)
         _fft_pointwise_sd_fft(sdctx, ii, jj, 0, trunc, limbs);
      else
#endif
      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         if (ii != jj) mpn_normmod_2expp1(jj[j], limbs);

         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
      }

      ifft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

      for (j = 0; j < trunc; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
         mpn_normmod_2expp1(ii[j], limbs);
      }
   } else
   {
      trunc = 2*sqrt*((trunc + 2*sqrt - 1)/(2*sqrt));

      /* The fused inner couples each pointwise batch to
         bandwidth-bound column transforms inside one work unit;
         measurement shows it scales well only when the unit count
         (2n/n1 columns plus tail rows) can feed the threads. Below
         that, separate full transforms with single-thread locality
         plus a flat threaded pointwise -- many small compute-heavy
         units, no barriers -- parallelize better. The factor of
         eight units per thread is a tuning parameter. */
#if FLINT_HAVE_FFT_SMALL
      if (sdctx != NULL
          && (2*n)/sqrt + (trunc - 2*n)/sqrt
             < 8*flint_get_num_threads())
      {
         fft_mfa_truncate_sqrt2(ii, n, w, t1, t2, s1, sqrt, trunc);

         if (ii != jj)
            fft_mfa_truncate_sqrt2(jj, n, w, t1, t2, s1, sqrt, trunc);

         _fft_pointwise_sd_fft(sdctx, ii, jj, 0, 2*n, limbs);
         _fft_pointwise_rows_sd_fft(sdctx, ii, jj, 2*n, sqrt,
                                    (trunc - 2*n)/sqrt,
                                    depth - depth/2 + 1, limbs);

         ifft_mfa_truncate_sqrt2(ii, n, w, t1, t2, s1, sqrt, trunc);

         /* the separate-transform path divides out the inverse
            scaling explicitly; the fused inner folds it into its
            workers */
         for (j = 0; j < trunc; j++)
         {
            mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
            mpn_normmod_2expp1(ii[j], limbs);
         }
      }
      else
#endif
      {
         fft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);

         if (ii != jj)
            fft_mfa_truncate_sqrt2_outer(jj, n, w, t1, t2, s1, sqrt,
                                         trunc);

#if FLINT_HAVE_FFT_SMALL
         if (sdctx != NULL)
            fft_mfa_truncate_sqrt2_inner_sd_fft(sdctx, ii, jj, n, w,
                                          t1, t2, s1, sqrt, trunc, tt);
         else
#endif
            fft_mfa_truncate_sqrt2_inner(ii, jj, n, w, t1, t2, s1,
                                         sqrt, trunc, tt);

         ifft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt,
                                       trunc);
      }
   }
}

/* Negacyclic-aware ring selection. The instantiated digit sizes
   form a lattice of rings N' = m b (m a power of two, b below),
   dense enough that some ring lies within a few percent of any
   target. Candidates are admitted when their ring growth over the
   standard Nussbaumer rounding stays within a per-prime-count cap
   -- low prime counts are robustly profitable and may pay more
   ring growth; high counts must come nearly free -- and the
   admitted candidate minimizing np * m * (depth + 4 + b/16)
   (transform work plus per-digit input/crt work) is selected.
   The caps and the cost weights are tuning parameters: they were
   calibrated conservatively against C-build kernel measurements
   and should be re-fit on tuned-assembly targets, where measured
   negmul advantages are larger and admit looser caps.
   On success *m_out receives the digit count and the returned limb
   count defines N' = limbs * FLINT_BITS; otherwise *m_out = 0 and
   the standard rounding is returned, keeping fft_mulmod_2expp1. */
slong fft_adjust_limbs_sd_fft(slong tight, slong n_outer, slong * m_out)
{
#if !FLINT_HAVE_FFT_SMALL
    *m_out = 0;
    return fft_adjust_limbs(tight);
#else

    static const slong bs[] = { 64, 84, 88, 92, 112, 116, 120, 126,
                                128, 136, 140, 144, 160, 164, 168,
                                184, 188, 192 };
    /* growth caps, indexed by np - 3 */
    static const double cap[6] = { 1.35, 1.25, 1.15, 1.12, 1.05, 1.02 };
    slong std = fft_adjust_limbs(tight);
    slong ii, bestl = std, bestm = 0;
    double bestc = 0.0;

    for (ii = 0; ii < (slong) (sizeof(bs) / sizeof(bs[0])); ii++)
    {
        slong b = bs[ii], depth, np;
        slong m = 16;
        double cost, grow;

        while (m * b < FLINT_BITS * tight)
            m <<= 1;
        if ((m * b) % FLINT_BITS != 0)
            continue;
        if (m * b < FLINT_FFT_SD_FFT_MIN_RING)
            continue;
        /* outer transform: the ring must be a multiple of the
           transform length with an even bits-per-point w; the digit
           count m is otherwise unrelated to n_outer */
        if ((m * b) % n_outer != 0 || ((m * b) / n_outer) % 2 != 0)
            continue;
        depth = FLINT_BIT_COUNT((ulong) m) - 1;
        np = sd_fft_mpn_mulmod_2expp1_choose_np(get_default_mpn_ctx(), b, depth);
        if (np < 0)
            continue;
        grow = (double) (m * b) / (double) (FLINT_BITS * std);
        if (grow > cap[np - 3])
            continue;
        cost = (double) np * (double) m
               * ((double) depth + 4.0 + (double) b / 16.0);
        if (bestm == 0 || cost < bestc)
        {
            bestc = cost;
            bestl = (m * b) / FLINT_BITS;
            bestm = m;
        }
    }
    *m_out = bestm;
    return bestl;
#endif
}
