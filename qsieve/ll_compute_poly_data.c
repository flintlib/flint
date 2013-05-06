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
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

void balance4(qs_t qs_inf, mp_limb_t * A_ind, 
     prime_t * factor_base, long min, long fact, long span, mp_limb_t target)
{
    long i, j;
    
    mp_limb_t prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p;
    i = fact;
    j = i + 1;
    while (j < min + span)
    {
        while (j < min + span)
        {
            if (prod*factor_base[i].p*factor_base[j].p >= target/P_GOODNESS)
                break;
            j++;
        }
        i++;
        j = i + 1;
    }
    A_ind[2] = i;
    A_ind[3] = j;
}

void balance5(qs_t qs_inf, mp_limb_t * A_ind, 
     prime_t * factor_base, long min, long high, long span, mp_limb_t target)
{
    long i, j;
    
    mp_limb_t prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p 
                   * factor_base[A_ind[2]].p;
    i = A_ind[2] + 1;
    j = i + 1;
    while (j < min + span)
    {
        while (j < min + span)
        {
            if (prod*factor_base[i].p*factor_base[j].p >= target/P_GOODNESS)
                break;
            j++;
        }
        i++;
        j = i + 1;
    }   
    A_ind[3] = i;
    A_ind[4] = j;
}

void try_compute_A(qs_t qs_inf)
{
   long min = qs_inf->min;
   long span = qs_inf->span;
   long fact = qs_inf->fact;
   long mid = qs_inf->mid;
   long high = qs_inf->high;
   long s = qs_inf->s;
   mp_limb_t * A_ind = qs_inf->A_ind;
   mp_limb_t target = qs_inf->target_A;
   prime_t * factor_base = qs_inf->factor_base;
   long i, j;
   mp_limb_t prod;
       
   if (qs_inf->A == 0) /* this is our first poly */
   {
       A_ind[0] = min;

       /* try to pick prime factors of A whose product is not much smaller than target_A */
       switch (s) /* we can only have up to 5 factors in A for small factorisations */
       {
       case 1:
           break;
       case 2:
           prod = factor_base[A_ind[0]].p;
           i = A_ind[0] + 1;
           while (prod*factor_base[i].p < target/P_GOODNESS2 && i + 1 < min + span) i++;
           A_ind[1] = i;
           break;
       case 3:
           A_ind[1] = mid;
           prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p;
           i = A_ind[1] + 1;
           while (prod*factor_base[i].p < target/P_GOODNESS && i + 1 < min + span) i++;
           A_ind[2] = i;
           break;
       case 4:
           A_ind[1] = A_ind[0] + 1;
           balance4(qs_inf, A_ind, factor_base, min, fact, span, target);
           break;
       case 5:
           A_ind[1] = A_ind[0] + 1;
           A_ind[2] = mid;
           balance5(qs_inf, A_ind, factor_base, min, high, span, target);
           break;
       }
   } else /* update to the next poly */
   {
       switch (s) 
       {
       case 1:
           if (A_ind[0] + 1 < min + span)
               A_ind[0]++;
           else
               goto out_of_polys;
           break;
       case 2:
           i = A_ind[0];
           j = A_ind[1] + 1;
           if (j < min + span && factor_base[i].p * factor_base[j].p < P_GOODNESS2*target)
               A_ind[1] = j; /* we can just increment second index */
           else /* must increment first index */
           {
               i++;
               if (i < fact) /* find first appropriate second index */
               {
                  A_ind[0] = i;
                  prod = factor_base[A_ind[0]].p;
                  i = A_ind[0] + 1;
                  while (prod*factor_base[i].p < target/P_GOODNESS && i + 1 < min + span) i++;
                  A_ind[1] = i;
               } else 
                   goto out_of_polys;
           }
           break;
       case 3:
           prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p;
           j = A_ind[2] + 1;
           if (j < min + span && prod * factor_base[j].p < P_GOODNESS*target)
               A_ind[2] = j; /* increment third index */
           else 
           {
               A_ind[1]++; /* increment second index */
               i = A_ind[1] + 1;
               prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p;
               if (i < min + span && prod * factor_base[i].p < P_GOODNESS*target)
               {
                   while (prod*factor_base[i].p < target/P_GOODNESS && i + 1 < min + span) i++;
                   A_ind[2] = i; 
               } else /* must increment first index */
               {
                   A_ind[0]++;
                   if (A_ind[0] < mid)
                   {
                       A_ind[1] = mid;
                       prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p;
                       i = A_ind[1] + 1;
                       while (prod*factor_base[i].p < target/P_GOODNESS && i + 1 < min + span) i++;
                       A_ind[2] = i;
                   } else
                       goto out_of_polys;
               }
           }
           break;
       case 4:
           prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p * factor_base[A_ind[2]].p;
           j = A_ind[3] + 1;
           if (j < min + span && prod * factor_base[j].p < P_GOODNESS*target)
               A_ind[3] = j; /* increment fourth index */
           else 
           {
               A_ind[2]++; /* increment third index */
               i = A_ind[2] + 1;
               prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p * factor_base[A_ind[2]].p;
               if (i < min + span && prod * factor_base[i].p < P_GOODNESS*target)
               {
                   while (prod*factor_base[i].p < target/P_GOODNESS && i + 1 < min + span) i++;
                   A_ind[3] = i; 
               } else
               {
                   A_ind[1]++; /* increment second index */
                   if (A_ind[1] < fact)
                   {
                       balance4(qs_inf, A_ind, factor_base, min, fact, span, target);
                   } else
                   {
                       A_ind[0]++; /* increment first factor */
                       A_ind[1] = A_ind[0] + 1;
                       if (A_ind[1] < fact)
                       {
                           balance4(qs_inf, A_ind, factor_base, min, fact, span, target);
                       } else
                          goto out_of_polys;
                   }
               }
           }
           break;
       case 5:
           prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p * factor_base[A_ind[2]].p * factor_base[A_ind[3]].p;
           j = A_ind[4] + 1;
           if (j < min + span && prod * factor_base[j].p < P_GOODNESS*target)
               A_ind[4] = j; /* increment fifth index */
           else 
           {
               A_ind[3]++; /* increment fourth index */
               i = A_ind[3] + 1;
               prod = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p * factor_base[A_ind[2]].p * factor_base[A_ind[3]].p;
               if (i < min + span && prod * factor_base[i].p < P_GOODNESS*target)
               {
                   while (prod*factor_base[i].p < target/P_GOODNESS && i + 1 < min + span) i++;
                   A_ind[4] = i; 
               } else
               {
                   A_ind[2]++; /* increment third index */
                   if (A_ind[2] < high)
                   {
                       balance5(qs_inf, A_ind, factor_base, min, high, span, target);
                   } else
                   {
                       A_ind[1]++; /* increment second index */
                       if (A_ind[1] < high)
                       {
                           A_ind[2] = mid;
                           balance5(qs_inf, A_ind, factor_base, min, high, span, target);
                       } else
                       {
                           A_ind[0]++; /* increment first factor */
                           A_ind[1] = A_ind[0] + 1;
                           if (A_ind[1] < mid)
                           {
                               A_ind[2] = mid;
                               balance5(qs_inf, A_ind, factor_base, min, high, span, target);
                           } else
                               goto out_of_polys;
                       }
                   }
               }
           }
           break;
       }
   }

   qs_inf->A = 1;
   for (i = 0; i < s; i++)
       qs_inf->A *= factor_base[A_ind[i]].p;

   return;

out_of_polys:
   printf("Out of polynomials, s = %ld\n", qs_inf->s);
   abort();
}

void qsieve_ll_compute_A(qs_t qs_inf)
{
    long i;
    
    do
    {
        try_compute_A(qs_inf);
    } while (((qs_inf->A > P_GOODNESS * qs_inf->target_A 
          || qs_inf->A < qs_inf->target_A / P_GOODNESS) && qs_inf->s > 2)
          || (((qs_inf->A > P_GOODNESS2 * qs_inf->target_A 
          || qs_inf->A < qs_inf->target_A / P_GOODNESS2) && qs_inf->s == 2)));

#if QS_DEBUG > 1
   printf("A = %ld, target A = %ld\n", qs_inf->A, qs_inf->target_A);
#endif    
 
    for (i = 0; i < qs_inf->s; i++)
    {
        mp_limb_t p = qs_inf->factor_base[qs_inf->A_ind[i]].p;
        qs_inf->inv_p2[i] = n_preinvert_limb(p*p);
    }      
}

void qsieve_ll_compute_B_terms(qs_t qs_inf)
{
   long s = qs_inf->s;
   mp_limb_t * A_ind = qs_inf->A_ind;
   mp_limb_t * A_modp = qs_inf->A_modp;
   mp_limb_t * B_terms = qs_inf->B_terms;
   prime_t * factor_base = qs_inf->factor_base;
   mp_limb_t A = qs_inf->A;
   mp_limb_t B;
   mp_limb_t p, temp, temp2, pinv;
   long i;
   
   for (i = 0; i < s; i++)
   {
      p = factor_base[A_ind[i]].p;
      pinv = factor_base[A_ind[i]].pinv;
      temp = A/p; /* TODO: possibly use precomputed inverse here */ 
      A_modp[i] = (temp2 = n_mod2_preinv(temp, p, pinv));
      temp2 = n_invmod(temp2, p);
      temp2 = n_mulmod2_preinv(temp2, qs_inf->sqrts[A_ind[i]], p, pinv);
      if (temp2 > p/2) temp2 = p - temp2;
      B_terms[i] = temp*temp2;     
   }
   
   B = B_terms[0];
   for (i = 1; i < s; i++)
   {
      B += B_terms[i];
   }
   qs_inf->B = B;
}

void qsieve_ll_compute_off_adj(qs_t qs_inf)
{
   long num_primes = qs_inf->num_primes;
   mp_limb_t A = qs_inf->A;
   mp_limb_t B = qs_inf->B;
   mp_limb_t * A_inv = qs_inf->A_inv;
   mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
   mp_limb_t * B_terms = qs_inf->B_terms;
   mp_limb_t * soln1 = qs_inf->soln1;
   mp_limb_t * soln2 = qs_inf->soln2;
   int * sqrts = qs_inf->sqrts;
   prime_t * factor_base = qs_inf->factor_base;
   long s = qs_inf->s;
   mp_limb_t p, temp, pinv;
   long i, j;
   
   for (i = 2; i < num_primes; i++) /* skip k and 2 */
   {
      p = factor_base[i].p;
      pinv = factor_base[i].pinv;
      
      A_inv[i] = n_invmod(n_mod2_preinv(A, p, pinv), p);
             
      for (j = 0; j < s; j++)
      {
         temp = n_mod2_preinv(B_terms[j], p, pinv);
         temp = n_mulmod2_preinv(temp, A_inv[i], p, pinv);
         temp *= 2;
         if (temp >= p) temp -= p;
         A_inv2B[j][i] = temp;
      }
             
      temp = n_mod2_preinv(B, p, pinv);
      temp = sqrts[i] + p - temp;
      temp *= A_inv[i];
      temp += qs_inf->sieve_size/2;
      soln1[i] = n_mod2_preinv(temp, p, pinv); 
      temp = p - sqrts[i];
      if (temp == p) temp -= p;
      temp = n_mulmod2_preinv(temp, A_inv[i], p, pinv);
      temp *= 2;
      if (temp >= p) temp -= p;      
      soln2[i] = temp + soln1[i];
      if (soln2[i] >= p) soln2[i] -= p;
   }  
}

void qsieve_ll_compute_A_factor_offsets(qs_t qs_inf)
{
   long s = qs_inf->s;
   mp_limb_t * A_ind = qs_inf->A_ind;
   mp_limb_t * A_modp = qs_inf->A_modp;
   mp_limb_t * soln1 = qs_inf->soln1;
   mp_limb_t * soln2 = qs_inf->soln2;
   mp_limb_t p, D;
   mp_limb_t hi = qs_inf->hi;
   mp_limb_t lo = qs_inf->lo;
   mp_limb_t B = qs_inf->B;
   mp_limb_t temp, temp2, B_modp2, index, p2; 
   prime_t * factor_base = qs_inf->factor_base;
   mp_limb_t * inv_p2 = qs_inf->inv_p2;
   mp_limb_t pinv;
   long j;
   
   for (j = 0; j < s; j++)
   {
      index = A_ind[j];
      p = factor_base[index].p;
      p2 = p*p;
      pinv = factor_base[index].pinv;
      D = n_ll_mod_preinv(hi, lo, p*p, inv_p2[j]);    
      if ((mp_limb_signed_t) B < 0) 
      {
          B_modp2 = n_mod2_preinv(-B, p2, inv_p2[j]);
          B_modp2 = p2 - B_modp2;
          if (B_modp2 == p2) B_modp2 = 0;
      } else
          B_modp2 = n_mod2_preinv(B, p2, inv_p2[j]);
      temp = B_modp2*A_modp[j];
      temp = n_mod2_preinv(temp, p, pinv); 
      temp2 = n_invmod(temp, p);
      D -= (B_modp2*B_modp2);
      if ((mp_limb_signed_t) D < 0) 
          temp = -(-D/p); /* TODO consider using precomputed inverse */
      else 
          temp = (D/p); /* TODO consider using precomputed inverse */
      temp *= temp2;
      temp += qs_inf->sieve_size/2;
      if ((mp_limb_signed_t) temp < 0) 
      {
         temp = p - n_mod2_preinv(-temp, p, pinv);
         if (temp == p) temp = 0;
      }
      else temp = n_mod2_preinv(temp, p, pinv);
      soln1[index] = temp;
      soln2[index] = -1;
   }
}          

void qsieve_ll_compute_C(qs_t qs_inf)
{
   mp_limb_t A = qs_inf->A;
   mp_limb_t B = qs_inf->B;
   
   if ((mp_limb_signed_t) B < 0L) B = -B;
   fmpz_set_ui(qs_inf->C, B);
   fmpz_mul_ui(qs_inf->C, qs_inf->C, B);
   fmpz_sub(qs_inf->C, qs_inf->C, qs_inf->kn);
   fmpz_divexact_ui(qs_inf->C, qs_inf->C, A);
} 

void qsieve_ll_compute_poly_data(qs_t qs_inf)
{
   qsieve_ll_compute_A(qs_inf);
   qsieve_ll_compute_B_terms(qs_inf);
   qsieve_ll_compute_off_adj(qs_inf);
   qsieve_ll_compute_A_factor_offsets(qs_inf);
   qsieve_ll_compute_C(qs_inf);        
}

