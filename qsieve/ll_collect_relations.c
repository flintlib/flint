/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2006, 2011 William Hart

******************************************************************************/

#undef ulong /* avoid clash with stdlib */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long 

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

void qsieve_ll_do_sieving(qs_t qs_inf, char * sieve)
{
   len_t num_primes = qs_inf->num_primes;
   mp_limb_t * soln1 = qs_inf->soln1;
   mp_limb_t * soln2 = qs_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   mp_limb_t p;
   char * end = sieve + qs_inf->sieve_size;
   register char * pos1;
   register char * pos2;
   register char * bound;  
   len_t size;
   len_t diff;
   len_t pind;
   
   memset(sieve, 0, qs_inf->sieve_size + sizeof(ulong));
   *end = (char) 255;
   
   for (pind = qs_inf->small_primes; pind < num_primes; pind++) 
   {
      if (soln2[pind] == -1) continue; /* don't sieve with A factors */
      
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

len_t qsieve_ll_evaluate_candidate(qs_t qs_inf, len_t i, char * sieve)
{
   len_t bits, exp, extra_bits;
   mp_limb_t modp, prime;
   len_t num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   fac_t * factor = qs_inf->factor;
   mp_limb_t * soln1 = qs_inf->soln1;
   mp_limb_t * soln2 = qs_inf->soln2;
   len_t * small = qs_inf->small;
   mp_limb_t A = qs_inf->A;
   mp_limb_t B = qs_inf->B;
   mp_limb_t pinv;
   len_t num_factors = 0;
   len_t relations = 0;
   len_t j;
   
   fmpz_t X, Y, res, p;
   fmpz_init(X); 
   fmpz_init(Y); 
   fmpz_init(res); 
   fmpz_init(p); 
    
   fmpz_set_ui(X, i);
   fmpz_sub_ui(X, X, qs_inf->sieve_size/2); /* X */
     
#if (QS_DEBUG & 32)
   printf("i = "); fmpz_print(X); printf("\n");
#endif

   fmpz_mul_ui(Y, X, A);
   if ((mp_limb_signed_t) B < 0) 
   {
      fmpz_sub_ui(Y, Y, -B);  /* Y = AX + B */
      fmpz_sub_ui(res, Y, -B);  
   } else 
   {
      fmpz_add_ui(Y, Y, B);
      fmpz_add_ui(res, Y, B);
   }
   fmpz_mul(res, res, X);  
   fmpz_add(res, res, qs_inf->C); /* res = AX^2 + 2BX + C */
           
   bits = FLINT_ABS(fmpz_bits(res));
   bits -= BITS_ADJUST; 
   extra_bits = 0;
   
   fmpz_set_ui(p, 2); /* divide out by powers of 2 */
   exp = fmpz_remove(res, res, p);

#if (QS_DEBUG & 8)
   if (exp) printf("2^%ld ", exp);
#endif

   extra_bits += exp;
   small[1] = exp;
     
   if (factor_base[0].p != 1) /* divide out powers of the multiplier */
   {
      fmpz_set_ui(p, factor_base[0].p);
      exp = fmpz_remove(res, res, p);
      if (exp) extra_bits += exp*qs_inf->factor_base[0].size;
      small[0] = exp;

#if (QS_DEBUG & 8)
      if (exp) printf("%d^%ld ", factor_base[0].p, exp); 
#endif
   } else small[0] = 0;
     
   for (j = 2; j < qs_inf->small_primes; j++) /* pull out small primes */
   {
      prime = factor_base[j].p;
      pinv = factor_base[j].pinv;
      modp = n_mod2_preinv(i, prime, pinv);
      if ((modp == soln1[j]) || (modp == soln2[j]))
      {
         fmpz_set_ui(p, prime);
         exp = fmpz_remove(res, res, p);
         if (exp) extra_bits += qs_inf->factor_base[j].size;
         small[j] = exp;

#if (QS_DEBUG & 8)
         if (exp) 
         {
             fmpz_print(p);
             printf("^%ld ", exp); 
         }
#endif
      } else small[j] = 0;
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
         if (soln2[j] != -1)
         {
            if ((modp == soln1[j]) || (modp == soln2[j]))
            {
               fmpz_set_ui(p, prime);
               exp = fmpz_remove(res, res, p);
#if (QS_DEBUG & 8)
               if (exp) 
               {
                   fmpz_print(p);
                   printf("^%ld ", exp); 
               }
#endif
               if (exp) 
               {
                  extra_bits += qs_inf->factor_base[j].size;
                  factor[num_factors].ind = j;
                  factor[num_factors++].exp = exp; 
               }
            }
         } else
         {
            fmpz_set_ui(p, prime);
            exp = fmpz_remove(res, res, p);
            factor[num_factors].ind = j;
            factor[num_factors++].exp = exp + 1; 

#if (QS_DEBUG & 8)
            if (exp) 
            {
                fmpz_print(p);
                printf("^%ld ", exp); 
            }
#endif
         }    
      }

      if (fmpz_cmp_ui(res, 1) == 0 || fmpz_cmp_si(res, -1) == 0) /* We've found a relation */
      {
         mp_limb_t * A_ind = qs_inf->A_ind;
         len_t i;

         for (i = 0; i < qs_inf->s; i++) /* Commit any outstanding A factors */
         {
            if (A_ind[i] >= j)
            {
               factor[num_factors].ind = A_ind[i];
               factor[num_factors++].exp = 1; 
            }
         }

         qs_inf->num_factors = num_factors;
         relations += qsieve_ll_insert_relation(qs_inf, Y);  /* Insert the relation in the matrix */
        
         if (qs_inf->num_relations >= qs_inf->buffer_size)
         {
            printf("Error: too many duplicate relations!\n");
            printf("s = %ld, bits = %ld\n", qs_inf->s, qs_inf->bits);
            abort();
         }

         goto cleanup;
      }
   }

#if (QS_DEBUG & 8)
   printf("\n");
#endif
  
cleanup:
   fmpz_clear(X);
   fmpz_clear(Y);
   fmpz_clear(res);
   fmpz_clear(p);
      
   return relations;
}

len_t qsieve_ll_evaluate_sieve(qs_t qs_inf, char * sieve)
{
   len_t i = 0, j = 0;
   ulong * sieve2 = (ulong *) sieve;
   char bits = qs_inf->sieve_bits;
   len_t rels = 0;

#if (QS_DEBUG & 16)
   len_t stats_limit;
   for (i = 0; i < 256; i++)
       qs_inf->sieve_tally[i] = 0;
#endif

#if (QS_DEBUG & 4)
   printf("%ldX^2+2*%ldX+", qs_inf->A, qs_inf->B);
   fmpz_print(qs_inf->C); printf("\n");
#endif

   while (j < qs_inf->sieve_size/sizeof(ulong))
   {
#if FLINT64
       while ((sieve2[j] & 0xE0E0E0E0E0E0E0E0UL) == 0) 
#else
       while ((sieve2[j] & 0xE0E0E0E0UL) == 0) 
#endif
       {
#if (QS_DEBUG & 16)
           for (i = j*sizeof(ulong); i < (j+1)*sizeof(ulong) && i < qs_inf->sieve_size; i++)
               qs_inf->sieve_tally[(int)sieve[i]]++;
#endif
           j++;
       }

       i = j*sizeof(ulong);

       while (i < (j+1)*sizeof(ulong) && i < qs_inf->sieve_size)
       {
#if (QS_DEBUG & 16)
           qs_inf->sieve_tally[(int)sieve[i]]++;
#endif
           if (sieve[i] > bits) 
               rels += qsieve_ll_evaluate_candidate(qs_inf, i, sieve);

           i++;
       }
       j++;
   }

#if (QS_DEBUG & 16)
   for (stats_limit = 255; stats_limit >= 0; stats_limit--)
       if (qs_inf->sieve_tally[stats_limit] != 0)
           break;
   
   for (i = 0; i <= stats_limit; i++)
   {
       if ((i % 16) == 0)
           printf("|%ld:", i);
       printf(" %ld", qs_inf->sieve_tally[i]);
   }
   printf("|\n");
   printf("Total of %ld relations for this sieve interval\n", rels);
#endif
 
   return rels;
}

void qsieve_ll_update_offsets(int poly_add, mp_limb_t * poly_corr, qs_t qs_inf)
{
   len_t num_primes = qs_inf->num_primes;
   mp_limb_t * soln1 = qs_inf->soln1;
   mp_limb_t * soln2 = qs_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   mp_limb_t p, correction;
   len_t pind;

   for (pind = 2; pind < num_primes; pind++) 
   {
      p = factor_base[pind].p;
      correction = (poly_add ? p - poly_corr[pind] : poly_corr[pind]);
      soln1[pind] += correction;
      if (soln1[pind] >= p) soln1[pind] -= p;
      if (soln2[pind] == -1) continue;
      soln2[pind] += correction;
      if (soln2[pind] >= p) soln2[pind] -= p; 
   }
}  

len_t qsieve_ll_collect_relations(qs_t qs_inf, char * sieve)
{
   len_t s = qs_inf->s;
   mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
   
   mp_limb_t * poly_corr;
   len_t relations = 0;
   len_t poly_index, j;
   int poly_add;
   
   qsieve_ll_compute_poly_data(qs_inf);
   
   for (poly_index = 1; poly_index < (1<<(s - 1)); poly_index++)
   {
      for (j = 0; j < s; j++)
         if (((poly_index >> j) & 1UL) != 0UL) break;
      
      poly_add = ((poly_index >> j) & 2);
      
      poly_corr = A_inv2B[j];
           
      qsieve_ll_do_sieving(qs_inf, sieve);
      
      relations += qsieve_ll_evaluate_sieve(qs_inf, sieve);
      
      qsieve_ll_update_offsets(poly_add, poly_corr, qs_inf);
      
      if (poly_add) qs_inf->B += (2*qs_inf->B_terms[j]); 
      else qs_inf->B -= (2*qs_inf->B_terms[j]); 
      
      qsieve_ll_compute_C(qs_inf);          
      
      qsieve_ll_compute_A_factor_offsets(qs_inf);    

      if (qs_inf->columns >= qs_inf->num_primes + qs_inf->extra_rels)
          break;
   }
   
   if (qs_inf->columns < qs_inf->num_primes + qs_inf->extra_rels)
   {
       qsieve_ll_do_sieving(qs_inf, sieve);
      
       relations += qsieve_ll_evaluate_sieve(qs_inf, sieve);

       relations += qsieve_ll_merge_relations(qs_inf);
   }
   
   return relations;
}
