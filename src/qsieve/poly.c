/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

/*
   Set up the data for sieving and checking candidates that depends only on
   the factor base: the ranges of primes sieved in the different ways (see the
   comment at BLOCK_SIZE in qsieve.h), the slices of bucket primes, the large
   prime bound, and the factor base data for trial division.
*/
static void
qsieve_sieve_init(qs_t qs_inf)
{
   slong num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   slong i, k, bp, lp;

   qs_inf->num_blocks = (qs_inf->sieve_size + QS_L1_SIZE - 1) / QS_L1_SIZE;

   for (lp = qs_inf->small_primes; lp < num_primes; lp++)
   {
      if (factor_base[lp].p >= QS_SMALL_BOUND)
         break;
   }

   for (bp = lp; bp < num_primes; bp++)
   {
      if (factor_base[bp].p >= BLOCK_SIZE)
         break;
   }

   qs_inf->l1_prime = lp;
   qs_inf->bucket_prime = bp;

   /* first prime that is at least the sieve size (the possible factors of A
      are much smaller, but make sure they are not included) */
   for (i = FLINT_MAX(bp, qs_inf->high); i < num_primes; i++)
   {
      if (factor_base[i].p >= (ulong) qs_inf->sieve_size)
         break;
   }

   qs_inf->huge_prime = i;

   qs_inf->lp_bound = FLINT_MIN(QS_LP_MULT * factor_base[num_primes - 1].p,
                                UWORD(0xffffffff));

   /* data for trial division of candidates, see qsieve_evaluate_candidate */
   qs_inf->fb_p32 = flint_malloc(num_primes*sizeof(uint32_t));
   qs_inf->fb_inv32 = flint_malloc(num_primes*sizeof(uint32_t));
   qs_inf->fb_lim32 = flint_malloc(num_primes*sizeof(uint32_t));

   for (i = 0; i < num_primes; i++)
   {
      ulong p = factor_base[i].p;

      /* the first three entries are the multiplier, 2 and -1 */
      if (i >= 3)
      {
         FLINT_ASSERT(p % 2 == 1 && p < UWORD(1) << 31);
         qs_inf->fb_p32[i] = p;
         qs_inf->fb_inv32[i] = (uint32_t) n_binvert(p);
         qs_inf->fb_lim32[i] = UWORD(0xffffffff) / p;
      }
      else
      {
         qs_inf->fb_p32[i] = qs_inf->fb_inv32[i] = qs_inf->fb_lim32[i] = 0;
      }
   }

   /* the bucket sieve is used unless the sieve is small */
   qs_inf->use_buckets = (qs_inf->sieve_size >= 2*QS_L1_SIZE);

   qs_inf->num_slices = 0;
   qs_inf->slice_start = NULL;
   qs_inf->slice_size = NULL;

   if (!qs_inf->use_buckets)
      return;

   /*
      A slice is a run of consecutive bucket primes of the same size, of
      length less than 2^16 so that the index within the slice fits in the
      top half of a bucket entry.
   */
   k = 0;
   qs_inf->slice_start = flint_malloc(sizeof(slong));
   for (i = bp; i < num_primes; i++)
   {
      if (i == bp || factor_base[i].size != factor_base[i - 1].size
                  || i - qs_inf->slice_start[k - 1] >= 65535)
      {
         qs_inf->slice_start = flint_realloc(qs_inf->slice_start, (k + 2)*sizeof(slong));
         qs_inf->slice_size = flint_realloc(qs_inf->slice_size, (k + 1)*sizeof(unsigned char));
         qs_inf->slice_start[k] = i;
         qs_inf->slice_size[k] = factor_base[i].size;
         k++;
      }
   }

   qs_inf->slice_start[k] = num_primes;
   qs_inf->num_slices = k;
}

/* allocate the buckets for one thread */
static void
qsieve_bucket_alloc(qs_t qs_inf, qs_poly_s * poly)
{
   slong b, i, alloc, num_blocks = qs_inf->num_blocks;
   double expected = 0.0;

   poly->bucket = NULL;
   poly->bucket_ptr = NULL;
   poly->bucket_alloc = NULL;
   poly->slice_end = NULL;

   if (!qs_inf->use_buckets)
      return;

   /* expected number of hits per block, plus some room for fluctuation */
   for (i = qs_inf->bucket_prime; i < qs_inf->num_primes; i++)
      expected += 2.0 * QS_L1_SIZE / qs_inf->factor_base[i].p;

   alloc = (slong) (1.1 * expected) + 4096;

   /* two more buckets, which receive the hits outside the sieve of the
      primes larger than it, see qsieve_do_sieving3 */
   poly->bucket = flint_malloc((num_blocks + 2)*sizeof(uint32_t *));
   poly->bucket_ptr = flint_malloc((num_blocks + 2)*sizeof(uint32_t *));
   poly->bucket_alloc = flint_malloc((num_blocks + 2)*sizeof(slong));
   poly->slice_end = flint_malloc((qs_inf->num_slices + 1)*num_blocks*sizeof(slong));

   for (b = 0; b < num_blocks + 2; b++)
   {
      poly->bucket[b] = flint_malloc(alloc*sizeof(uint32_t));
      poly->bucket_ptr[b] = poly->bucket[b];
      poly->bucket_alloc[b] = alloc;
   }
}

static void
qsieve_bucket_free(qs_t qs_inf, qs_poly_s * poly)
{
   slong b;

   if (poly->bucket == NULL)
      return;

   for (b = 0; b < qs_inf->num_blocks + 2; b++)
      flint_free(poly->bucket[b]);

   flint_free(poly->bucket);
   flint_free(poly->bucket_ptr);
   flint_free(poly->bucket_alloc);
   flint_free(poly->slice_end);

   poly->bucket = NULL;
}

ulong qsieve_poly_init(qs_t qs_inf)
{
   ulong num_primes = qs_inf->num_primes;
   ulong s = qs_inf->s;    /* number of prime factors in A coeff */
   uint32_t ** A_inv2B;
   slong i;

   fmpz_init(qs_inf->A);
   fmpz_init(qs_inf->B);
   fmpz_init(qs_inf->upp_bound);
   fmpz_init(qs_inf->low_bound);

   qs_inf->curr_subset = flint_malloc(s * sizeof(ulong));
   qs_inf->first_subset = flint_malloc(s * sizeof(ulong));
   qs_inf->B_terms = flint_malloc(s * sizeof(ulong));
   qs_inf->A_ind = flint_malloc(s * sizeof(ulong));
   qs_inf->A_divp = flint_malloc(s * sizeof(ulong));
   qs_inf->B0_terms = flint_malloc(s * sizeof(ulong));

   qs_inf->A_inv2B = flint_malloc(s * sizeof(uint32_t *));

   qs_inf->A_inv = flint_malloc(num_primes * sizeof(ulong));
   qs_inf->soln1 = flint_malloc(num_primes * sizeof(ulong));
   qs_inf->soln2 = flint_malloc(num_primes * sizeof(ulong));

   qs_inf->poly = flint_malloc((qs_inf->num_handles + 1)* sizeof(qs_poly_s));

   qsieve_sieve_init(qs_inf);

   for (i = 0; i <= qs_inf->num_handles ; i++)
   {
      fmpz_init(qs_inf->poly[i].B);
      fmpz_init(qs_inf->poly[i].C);
      qs_inf->poly[i].C_valid = 0;
      qs_inf->poly[i].A_root = flint_malloc(s*sizeof(ulong));
      qs_inf->poly[i].posn1 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].posn2 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].soln1 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].soln2 = flint_malloc((num_primes + 16)*sizeof(ulong));
      qs_inf->poly[i].small = flint_malloc(qs_inf->small_primes*sizeof(ulong));
      qs_inf->poly[i].factor = flint_malloc(qs_inf->max_factors*sizeof(fac_t));
      qsieve_bucket_alloc(qs_inf, qs_inf->poly + i);
   }

   A_inv2B = qs_inf->A_inv2B;

   for (i = 0; i < s; i++)
      A_inv2B[i] = flint_malloc(num_primes * sizeof(uint32_t));

   for (i = 0; i < s; i++)
   {
       fmpz_init(qs_inf->A_divp[i]);
       fmpz_init(qs_inf->B_terms[i]);
   }

   return 0;
}

void qsieve_poly_clear(qs_t qs_inf)
{
   slong i;

   /* already cleared (e.g. a small factor was found right after the factor
      base was increased) */
   if (qs_inf->poly == NULL)
      return;

   fmpz_clear(qs_inf->A);
   fmpz_clear(qs_inf->B);
   fmpz_clear(qs_inf->target_A);
   fmpz_clear(qs_inf->low_bound);
   fmpz_clear(qs_inf->upp_bound);

   for (i = 0; i < qs_inf->s; i++)
   {
       fmpz_clear(qs_inf->A_divp[i]);
       fmpz_clear(qs_inf->B_terms[i]);
   }

   flint_free(qs_inf->B_terms);
   flint_free(qs_inf->A_ind);
   flint_free(qs_inf->A_divp);
   flint_free(qs_inf->B0_terms);
   flint_free(qs_inf->A_inv);
   flint_free(qs_inf->soln1);
   flint_free(qs_inf->soln2);
   flint_free(qs_inf->curr_subset);
   flint_free(qs_inf->first_subset);

   if (qs_inf->A_inv2B != NULL)
   {
       for (i = 0; i < qs_inf->s; i++)
       {
           flint_free(qs_inf->A_inv2B[i]);
       }
   }

   flint_free(qs_inf->A_inv2B);

   for (i = 0; i <= qs_inf->num_handles; i++)
   {
      fmpz_clear(qs_inf->poly[i].B);
      fmpz_clear(qs_inf->poly[i].C);
      flint_free(qs_inf->poly[i].A_root);
      flint_free(qs_inf->poly[i].posn1);
      flint_free(qs_inf->poly[i].posn2);
      flint_free(qs_inf->poly[i].soln1);
      flint_free(qs_inf->poly[i].soln2);
      flint_free(qs_inf->poly[i].small);
      flint_free(qs_inf->poly[i].factor);
      qsieve_bucket_free(qs_inf, qs_inf->poly + i);
   }
   flint_free(qs_inf->poly);
   qs_inf->poly = NULL;

   flint_free(qs_inf->slice_start);
   flint_free(qs_inf->slice_size);
   flint_free(qs_inf->fb_p32);
   flint_free(qs_inf->fb_inv32);
   flint_free(qs_inf->fb_lim32);
   qs_inf->fb_p32 = qs_inf->fb_inv32 = qs_inf->fb_lim32 = NULL;
   qs_inf->slice_start = NULL;
   qs_inf->slice_size = NULL;
   qs_inf->num_slices = 0;
   qs_inf->use_buckets = 0;

   qs_inf->B_terms = NULL;
   qs_inf->A_ind = NULL;
   qs_inf->A_divp = NULL;
   qs_inf->B0_terms = NULL;
   qs_inf->A_inv = NULL;
   qs_inf->soln1 = NULL;
   qs_inf->soln2 = NULL;
   qs_inf->A_inv2B = NULL;
   qs_inf->curr_subset = NULL;
}
