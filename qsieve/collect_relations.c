/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qsieve.h"

#include <time.h>

void qsieve_do_sieving_old(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
   unsigned long num_primes = qs_inf->num_primes;
   int * soln1 = poly->soln1;
   int * soln2 = poly->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   ulong sieve_size = qs_inf->sieve_size;
   unsigned char * end = sieve + sieve_size;
   ulong sieve_fill = qs_inf->sieve_fill;
   ulong small_primes = qs_inf->small_primes;
   unsigned char * bound;
   unsigned char * pos1;
   unsigned char * pos2;
   unsigned long size;
   ulong p;
   ulong prime;

   const ulong second_prime = qs_inf->second_prime;
   
   memset(sieve, sieve_fill, sieve_size);
   *end = 255;
   
   for (prime = small_primes; prime < second_prime; prime++) 
   {
      if (soln2[prime] == 0) continue;
      
      p = factor_base[prime].p;
      size = factor_base[prime].size;
      pos1 = sieve + soln1[prime];
      pos2 = sieve + soln2[prime];
      bound = end - p;
        
      while (bound - pos1 > 0)  
      {  
         (*pos1)+=size, pos1+=p, (*pos2)+=size, pos2+=p;
      }
      if ((end - pos1 > 0) && (end - pos2 > 0))
      { 
         (*pos1)+=size, pos1+=p, (*pos2)+=size, pos2+=p;
      }
      if (end - pos2 > 0)
      { 
         (*pos2)+=size;
      }
      if (end - pos1 > 0)
      { 
         (*pos1)+=size;
      }
   }
   
   for (prime = second_prime; prime < num_primes; prime++) 
   {
      p = factor_base[prime].p;
      size = factor_base[prime].size;
      pos1 = sieve + soln1[prime];
      pos2 = sieve + soln2[prime];
        
      if (end - pos2 > 0)
      { 
         (*pos2)+=size;
      }
      if (end - pos1 > 0)
      { 
         (*pos1)+=size;
      }
   }
}

void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
   slong num_primes = qs_inf->num_primes;
   int * soln1 = poly->soln1;
   int * soln2 = poly->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   mp_limb_t p;

   unsigned char * end = sieve + qs_inf->sieve_size;
   register unsigned char * pos1;
   register unsigned char * pos2;
   register unsigned char * bound;
   slong size;
   slong diff;
   slong pind;

   memset(sieve, qs_inf->sieve_fill, qs_inf->sieve_size + sizeof(ulong));
   *end = (char) 255;

   for (pind = qs_inf->small_primes; pind < num_primes; pind++)
   {
      if (soln2[pind] == 0) continue; /* don't sieve with A factors */

      p = factor_base[pind].p;
      size = factor_base[pind].size;
      pos1 = sieve + soln1[pind];
      pos2 = sieve + soln2[pind];
      diff = pos2 - pos1;
      bound = end - 2*p;

      while (bound - pos1 > 0)
      {
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
      }

      while ((end - pos1 > 0) && (end - pos1 - diff > 0))
      {
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
      }

      pos2 = pos1 + diff;

      if (end - pos2 > 0)
      {
         (*pos2) += size;
      }

      if (end - pos1 > 0)
      {
         (*pos1) += size;
      }
   }
}

/*
   sieving routine, breaks sieve array into blocks then sieve
   each block, for whole factor base at once
*/

void qsieve_do_sieving2(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
    slong b, d1, d2, i;
    slong pind, size;
    mp_limb_t p;
    slong num_primes = qs_inf->num_primes;
    int * soln1 = poly->soln1;
    int * soln2 = poly->soln2;
    int * posn1 = poly->posn1;
    int * posn2 = poly->posn2;
    prime_t * factor_base = qs_inf->factor_base;

    unsigned char * B;
    register unsigned char * Bp;
    register unsigned char * pos;

    memset(sieve, qs_inf->sieve_fill, qs_inf->sieve_size + sizeof(ulong));
    sieve[qs_inf->sieve_size] = (char) 255;

    for (i = 0; i < num_primes; i++)
    {
        posn1[i] = soln1[i];
        posn2[i] = soln2[i] - posn1[i];
    }

    for (b = 1; b <= qs_inf->sieve_size / BLOCK_SIZE; b++)
    {
        B = sieve + b * BLOCK_SIZE;

        for (pind = qs_inf->small_primes; pind < qs_inf->second_prime; pind++)
        {
            if (soln2[pind] == 0)
                continue;

            p = factor_base[pind].p;
            size = factor_base[pind].size;
            d1 = posn2[pind];
            d2 = p - d1;
            Bp = B - 2*d1 - d2;
            pos = sieve + posn1[pind];

            while (pos < Bp)
            {
                (*pos) += size, (*(pos + d1)) += size, pos += p;
                (*pos) += size, (*(pos + d1)) += size, pos += p;
            }

            Bp = B - d1;

            while (pos < Bp)
            {
                (*pos) += size, 
                (*(pos + d1)) += size, pos += p;
            }

            if (pos < B)
            {
                (*pos) += size, pos += d1;
                posn2[pind] = d2;
            }
            else { posn2[pind] = d1; }

            posn1[pind] = (pos - sieve);
        }

        for (pind = qs_inf->second_prime; pind < num_primes; pind++)
        {
            p = factor_base[pind].p;

            if (soln2[pind] == 0)
                continue;

            size = factor_base[pind].size;
            pos = sieve + posn1[pind];

            if (pos < B)
            {
                (*pos) += size;
                pos += posn2[pind];

                if (pos < B)
                {
                    (*pos) += size;
                    pos += p - posn2[pind];
                }
                else { posn2[pind] = p - posn2[pind]; }

                posn1[pind] = (pos - sieve);
            }
            else { posn1[pind] = (pos - sieve); }
        }
    }
}

/* check position 'i' in sieve array for smoothness */

slong qsieve_evaluate_candidate(qs_t qs_inf, ulong i, unsigned char * sieve, qs_poly_t poly)
{
   slong bits, exp, extra_bits;
   mp_limb_t modp, prime;
   slong num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   slong * small = poly->small;
   fac_t * factor = poly->factor;
   int * soln1 = poly->soln1;
   int * soln2 = poly->soln2;
   mp_limb_t * A_ind = qs_inf->A_ind;
   mp_limb_t pinv;
   slong num_factors = 0;
   slong relations = 0;
   slong j, k;

   fmpz_t X, Y, C, res, p;

   fmpz_init(X);
   fmpz_init(Y);
   fmpz_init(res);
   fmpz_init(p);
   fmpz_init(C);

   qsieve_compute_C(C, qs_inf, poly);   
      
   fmpz_set_si(X, i - qs_inf->sieve_size / 2); /* X */

   fmpz_mul(Y, X, qs_inf->A);
   fmpz_add(Y, Y, poly->B); /* Y = AX+B */
   fmpz_add(res, Y, poly->B); /* Y = AX+2B */

   fmpz_mul(res, res, X);
   fmpz_add(res, res, C); /* res = AX^2 + 2BX + C */

#if QS_DEBUG & 128
   printf("res = "); fmpz_print(res); printf("\n");
   flint_printf("Poly: "); fmpz_print(qs_inf->A); flint_printf("*x^2 + 2*");
   fmpz_print(poly->B); flint_printf("*x + "); fmpz_print(C); printf("\n");
   flint_printf("x = %wd\n", i - qs_inf->sieve_size / 2);
#endif

   sieve[i] -= qs_inf->sieve_fill;
   bits = FLINT_ABS(fmpz_bits(res));
   bits -= BITS_ADJUST;
   extra_bits = 0;

   if (factor_base[0].p != 1) /* divide out powers of the multiplier */
   {
      fmpz_set_ui(p, factor_base[0].p);
      exp = fmpz_remove(res, res, p);
      if (exp) extra_bits += exp*qs_inf->factor_base[0].size;
      small[0] = exp;
#if QS_DEBUG & 128
      if (exp)
          flint_printf("%ld^%ld ", factor_base[0].p, exp);
#endif
   }
   else small[0] = 0;

   fmpz_set_ui(p, 2); /* divide out by powers of 2 */
   exp = fmpz_remove(res, res, p);
#if QS_DEBUG & 128
   if (exp)
       flint_printf("%ld^%ld ", 2, exp);
#endif

   extra_bits += exp;
   small[1] = exp;

   for (j = 3; j < qs_inf->small_primes; j++) /* pull out small primes */
   {
      prime = factor_base[j].p;
      pinv = factor_base[j].pinv;
      modp = n_mod2_preinv(i, prime, pinv);

      if (modp == soln1[j] || modp == soln2[j])
      {
         fmpz_set_ui(p, prime);
         exp = fmpz_remove(res, res, p);
         if (exp) extra_bits += qs_inf->factor_base[j].size;
         small[j] = exp;
#if QS_DEBUG & 128
         if (exp)
            flint_printf("%ld^%ld ", prime, exp);
#endif
      }
      else small[j] = 0;
   }

   if (extra_bits + sieve[i] > bits)
   {
      sieve[i] += extra_bits;

      /* pull out remaining primes */
      for (j = qs_inf->small_primes; j < num_primes && extra_bits < sieve[i]; j++)
      {
         prime = factor_base[j].p;
         pinv = factor_base[j].pinv;
         modp = n_mod2_preinv(i, prime, pinv);

         if (soln2[j] != 0)
         {
            if (modp == soln1[j] || modp == soln2[j])
            {
               fmpz_set_ui(p, prime);
               exp = fmpz_remove(res, res, p);
               extra_bits += qs_inf->factor_base[j].size*exp;
               factor[num_factors].ind = j;
               factor[num_factors++].exp = exp;
#if QS_DEBUG & 128
                  flint_printf("%ld^%ld ", prime, exp);
#endif

            }
         }
         else
         {
            fmpz_set_ui(p, prime);
            exp = fmpz_remove(res, res, p);
            factor[num_factors].ind = j;
            factor[num_factors++].exp = exp + 1;
#if QS_DEBUG & 128
            if (exp)
               flint_printf("%ld^%ld ", prime, exp);
#endif

         }
      }

#if QS_DEBUG & 128
      if (num_factors > 0)
         flint_printf("\n");
#endif
      if (fmpz_cmp_ui(res, 1) == 0 || fmpz_cmp_si(res, -1) == 0) /* We've found a relation */
      {
#if QS_DEBUG
         if (qs_inf->full_relation % 100 == 0)
            printf("%ld relations\n", qs_inf->full_relation);
#endif
         
         if (fmpz_cmp_si(res, -1) == 0)
            small[2] = 1;
         else
            small[2] = 0;

         for (k = 0; k < qs_inf->s; k++) /* Commit any outstanding A factors */
         {
            if (A_ind[k] >= j)
            {
               factor[num_factors].ind = A_ind[k];
               factor[num_factors++].exp = 1;
            }
         }
         
         factor[num_factors].ind = qs_inf->q_idx;
         factor[num_factors++].exp = 1;

         poly->num_factors = num_factors;

#pragma omp critical
         {
            qsieve_write_to_file(qs_inf, 1, Y, poly);
         
            qs_inf->full_relation++;
         }
         relations++;

#if 0
         if (small[2] != 0)
            flint_printf("-");

         if (small[0] != 0)
            flint_printf("%wd^%d * ", qs_inf->k, small[0]);

         if (small[1] != 0)
            flint_printf("%wd^%d * ", WORD(2), small[1]);

         for (k = 3; k < qs_inf->small_primes; k++)
         {
            if (small[k] != 0)
               flint_printf("%wd^%d * ", factor_base[k].p, small[k]);
         }

         for (k = 0; k < num_factors - 1; k++)
            flint_printf("%wd^%d * ", factor_base[factor[k].ind].p, factor[k].exp);

         if (num_factors > 0)
            flint_printf("%wd^%d\n", factor_base[factor[k].ind].p, factor[k].exp);
         else
            flint_printf("\n");
#endif
      }
      else
      {
          if (fmpz_sgn(res) < 0)
          {
              fmpz_neg(res, res);
              small[2] = 1;
          } else
              small[2] = 0;

          if (fmpz_bits(res) <= 30)
          {
              prime = fmpz_get_ui(res);
              
              if (prime < 60 * factor_base[qs_inf->num_primes - 1].p && n_gcd(prime, qs_inf->k) == 1)
              {
                  for (k = 0; k < qs_inf->s; k++)  /* commit any outstanding A factor */
                  {
                      if (A_ind[k] >= j)
                      {
                          factor[num_factors].ind = A_ind[k];
                          factor[num_factors++].exp = 1;
                      }
                  }

                  factor[num_factors].ind = qs_inf->q_idx;
                  factor[num_factors++].exp = 1;

                  poly->num_factors = num_factors;

#pragma omp critical
                  {
                     /* store this partial in file */

                     qsieve_write_to_file(qs_inf, prime, Y, poly);

                     qs_inf->edges++;

                     qsieve_add_to_hashtable(qs_inf, prime);
                  }
              }
          }
      }

      goto cleanup;
   }

cleanup:
   fmpz_clear(X);
   fmpz_clear(Y);
   fmpz_clear(C);
   fmpz_clear(res);
   fmpz_clear(p);

   return relations;
}

/* scan sieve array for possible candidate for smooth relations */

slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
    slong i = 0, j = 0;
    ulong * sieve2 = (ulong *) sieve;
    unsigned char bits = qs_inf->sieve_bits;
    slong rels = 0;

    while (j < qs_inf->sieve_size / sizeof(ulong))
    {
#if FLINT64
        while ((sieve2[j] & UWORD(0xC0C0C0C0C0C0C0C0)) == 0)
#else
        while ((sieve2[j] & UWORD(0xC0C0C0C0)) == 0)
#endif
        {
            j++;
        }
        
        i = j * sizeof(ulong);

        while (i < (j + 1) * sizeof(ulong) && i < qs_inf->sieve_size)
        {
            if (sieve[i] > bits)
               rels += qsieve_evaluate_candidate(qs_inf, i, sieve, poly);

            i++;
        }
        
        j++;
    }
    
    return rels;
}

/* procedure to call polynomial initialization and sieving procedure */

slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve)
{
    slong relations = 0, rels, i, j = 0;

    qsieve_init_poly_first(qs_inf);

#pragma omp parallel for
    for (i = 0; i < (1 << qs_inf->s); i++)
    {
#if HAVE_OPENMP
        unsigned char * thread_sieve = sieve + (qs_inf->sieve_size + sizeof(ulong) + 64)*omp_get_thread_num();
        qs_poly_s * thread_poly = qs_inf->poly + omp_get_thread_num();
#else
        unsigned char * thread_sieve = sieve;
        qs_poly_s * thread_poly = qs_inf->poly;
#endif

#pragma omp critical
        {
           if (j == 0)
              qsieve_poly_copy(thread_poly, qs_inf);
           else
           {
              qsieve_init_poly_next(qs_inf, j);
              qsieve_poly_copy(thread_poly, qs_inf);
           }
           j++;
        }

        if (qs_inf->sieve_size < 2*BLOCK_SIZE)
           qsieve_do_sieving(qs_inf, thread_sieve, thread_poly);
        else
           qsieve_do_sieving2(qs_inf, thread_sieve, thread_poly);

        rels = qsieve_evaluate_sieve(qs_inf, thread_sieve, thread_poly);

#pragma omp atomic
        relations += rels;
    }

    return relations;
}
