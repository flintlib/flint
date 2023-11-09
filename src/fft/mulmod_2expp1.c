/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fft.h"
#include "ulong_extras.h"
#include "fft_tuning.h"
#include "mpn_extras.h"

static mp_size_t mulmod_2expp1_table_n[FFT_N_NUM] = MULMOD_TAB;

void fft_naive_convolution_1(mp_limb_t * r, mp_limb_t * ii, mp_limb_t * jj, mp_size_t m)
{
   mp_size_t i, j;

   for (i = 0; i < m; i++)
      r[i] = ii[0]*jj[i];

   for (i = 1; i < m; i++)
   {
      for (j = 0; j < m - i; j++)
         r[i+j] += ii[i]*jj[j];

      for ( ; j < m; j++)
         r[i+j-m] -=ii[i]*jj[j];
   }
}

void _fft_mulmod_2expp1(mp_limb_t * r1, mp_limb_t * i1, mp_limb_t * i2,
                 mp_size_t r_limbs, flint_bitcnt_t depth, flint_bitcnt_t w)
{
   mp_size_t n = (UWORD(1)<<depth);
   flint_bitcnt_t bits1 = (r_limbs*FLINT_BITS)/(2*n);

   mp_size_t limb_add, limbs = (n*w)/FLINT_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, ll;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *s1, *r, *ii0, *jj0;
   mp_limb_t c;

   ii = flint_malloc((2*(n + n*size) + 4*n + 5*size)*sizeof(mp_limb_t));
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size)
   {
      ii[i] = ptr;
   }
   ii0 = ptr;
   t1 = ii0 + 2*n;
   t2 = t1 + size;
   s1 = t2 + size;
   r = s1 + size;
   tt = r + 2*n;

   if (i1 != i2)
   {
      jj = flint_malloc((2*(n + n*size) + 2*n)*sizeof(mp_limb_t));
      for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size)
      {
         jj[i] = ptr;
      }
      jj0 = ptr;
   } else
   {
      jj = ii;
      jj0 = ii0;
   }

   j = fft_split_bits(ii, i1, r_limbs, bits1, limbs);
   for ( ; j < 2*n; j++)
      flint_mpn_zero(ii[j], limbs + 1);

   for (i = 0; i < 2*n; i++)
      ii0[i] = ii[i][0];

   fft_negacyclic(ii, n, w, &t1, &t2, &s1);
   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   if (i1 != i2)
   {
      j = fft_split_bits(jj, i2, r_limbs, bits1, limbs);
      for ( ; j < 2*n; j++)
         flint_mpn_zero(jj[j], limbs + 1);

      for (i = 0; i < 2*n; i++)
         jj0[i] = jj[i][0];

      fft_negacyclic(jj, n, w, &t1, &t2, &s1);
   }

   for (j = 0; j < 2*n; j++)
   {
      if (i1 != i2) mpn_normmod_2expp1(jj[j], limbs);
      c = 2*ii[j][limbs] + jj[j][limbs];
      ii[j][limbs] = flint_mpn_mulmod_2expp1_basecase(ii[j], ii[j], jj[j], c, n*w, tt);
   }

   ifft_negacyclic(ii, n, w, &t1, &t2, &s1);

   fft_naive_convolution_1(r, ii0, jj0, 2*n);

   for (j = 0; j < 2*n; j++)
   {
      mp_limb_t t, cy2;

      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 1);
      mpn_normmod_2expp1(ii[j], limbs);

      t = ii[j][limbs];
      ii[j][limbs] = r[j] - ii[j][0];
      cy2 = mpn_add_1(ii[j], ii[j], limbs + 1, ii[j][limbs]);
      add_ssaaaa(r[j], ii[j][limbs], 0, ii[j][limbs], 0, t);
      if (cy2) r[j]++;
   }

   flint_mpn_zero(r1, r_limbs + 1);
   fft_combine_bits(r1, ii, 2*n - 1, bits1, limbs + 1, r_limbs + 1);

   /*
      as the negacyclic convolution has effectively done subtractions
      some of the coefficients will be negative, so need to subtract p
   */
   ll = 0;
   limb_add = bits1/FLINT_BITS;

   for (j = 0; j < 2*n - 2; j++)
   {
      if (r[j])
         mpn_sub_1(r1 + ll + 1, r1 + ll + 1, r_limbs - ll, 1);
      else if ((mp_limb_signed_t) ii[j][limbs] < 0) /* coefficient was -ve */
      {
         mpn_sub_1(r1 + ll + 1, r1 + ll + 1, r_limbs - ll, 1);
         mpn_sub_1(r1 + ll + limbs + 1, r1 + ll + limbs + 1, r_limbs - limbs - ll, 1);
      }

      ll += limb_add;
   }
   /* penultimate coefficient, top bit was already ignored */
   if (r[j] || (mp_limb_signed_t) ii[j][limbs] < 0) /* coefficient was -ve */
      mpn_sub_1(r1 + ll + 1, r1 + ll + 1, r_limbs - ll, 1);

   /* final coefficient wraps around */
   if (limb_add)
      r1[r_limbs] += mpn_add_n(r1 + r_limbs - limb_add, r1 + r_limbs - limb_add, ii[2*n - 1], limb_add);
   c = mpn_sub_n(r1, r1, ii[2*n - 1] + limb_add, limbs + 1 - limb_add);
   mpn_addmod_2expp1_1(r1 + limbs + 1 - limb_add, r_limbs - limbs - 1 + limb_add, -c);
   mpn_normmod_2expp1(r1, r_limbs);

   flint_free(ii);
   if (i1 != i2) flint_free(jj);
}

void fft_mulmod_2expp1(mp_limb_t * r, mp_limb_t * i1, mp_limb_t * i2,
                           mp_size_t n, mp_size_t w, mp_limb_t * tt)
{
   mp_size_t bits = n*w;
   mp_size_t limbs = bits/FLINT_BITS;
   flint_bitcnt_t depth1, depth = 1;

   mp_size_t w1, off;

   mp_limb_t c = 2*i1[limbs] + i2[limbs];

   if (c & 1)
   {
      mpn_neg(r, i1, limbs + 1);
      mpn_normmod_2expp1(r, limbs);
      return;
   } else if (c & 2)
   {
      mpn_neg(r, i2, limbs + 1);
      mpn_normmod_2expp1(r, limbs);
      return;
   }

   if (limbs <= FFT_MULMOD_2EXPP1_CUTOFF)
   {
      r[limbs] = flint_mpn_mulmod_2expp1_basecase(r, i1, i2, c, bits, tt);
      return;
   }

   while ((UWORD(1)<<depth) < bits) depth++;

   if (depth < 12) off = mulmod_2expp1_table_n[0];
   else off = mulmod_2expp1_table_n[FLINT_MIN(depth, FFT_N_NUM + 11) - 12];
   depth1 = depth/2 - off;

   w1 = bits/(UWORD(1)<<(2*depth1));

   _fft_mulmod_2expp1(r, i1, i2, limbs, depth1, w1);
}

slong fft_adjust_limbs(mp_size_t limbs)
{
   mp_size_t bits1 = limbs*FLINT_BITS, bits2;
   mp_size_t depth = 1, limbs2, depth1 = 1, depth2 = 1, adj;
   mp_size_t off1, off2;

   if (limbs <= FFT_MULMOD_2EXPP1_CUTOFF) return limbs;

   depth = FLINT_CLOG2(limbs);
   limbs2 = (WORD(1)<<depth); /* within a factor of 2 of limbs */
   bits2 = limbs2*FLINT_BITS;

   depth1 = FLINT_CLOG2(bits1);
   if (depth1 < 12) off1 = mulmod_2expp1_table_n[0];
   else off1 = mulmod_2expp1_table_n[FLINT_MIN(depth1, FFT_N_NUM + 11) - 12];
   depth1 = depth1/2 - off1;

   depth2 = FLINT_CLOG2(bits2);
   if (depth2 < 12) off2 = mulmod_2expp1_table_n[0];
   else off2 = mulmod_2expp1_table_n[FLINT_MIN(depth2, FFT_N_NUM + 11) - 12];
   depth2 = depth2/2 - off2;

   depth1 = FLINT_MAX(depth1, depth2);
   adj = (WORD(1)<<(depth1 + 1));
   limbs2 = adj*((limbs + adj - 1)/adj); /* round up number of limbs */
   bits1 = limbs2*FLINT_BITS;
   bits2 = (WORD(1)<<(depth1*2));
   bits1 = bits2*((bits1 + bits2 - 1)/bits2); /* round up bits */
   limbs = bits1/FLINT_BITS;

   return limbs;
}
