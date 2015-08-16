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

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

#include <time.h>

/* taken from FLINT2 */

void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve)
{
   slong num_primes = qs_inf->num_primes;
   mp_limb_t * soln1 = qs_inf->soln1;
   mp_limb_t * soln2 = qs_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   mp_limb_t p;
   char * end = sieve + qs_inf->sieve_size;
   register char * pos1;
   register char * pos2;
   register char * bound;
   slong size;
   slong diff;
   slong pind;

   memset(sieve, 10, qs_inf->sieve_size + sizeof(ulong));
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

void qsieve_do_sieving2(qs_t qs_inf, unsigned char * sieve)
{
    slong b, d1, d2, i, j;
    slong pind, size, r0;
    slong num_primes = qs_inf->num_primes;
    mp_limb_t p;
    mp_limb_t * soln1 = qs_inf->soln1;
    mp_limb_t * soln2 = qs_inf->soln2;
    mp_limb_t * xr1 = qs_inf->xr1;
    mp_limb_t * xr2 = qs_inf->xr2;
    prime_t * factor_base = qs_inf->factor_base;

    unsigned char * B;
    unsigned char * Bp;
    register unsigned char * x;

    memset(sieve, 0, qs_inf->sieve_size + sizeof(ulong));

    for (i = 0; i < num_primes; i++)
    {
        xr1[i] = soln1[i];
        xr2[i] = soln2[i] - xr1[i];
    }

    for (i = 0; i < num_primes; i++)
    {
        if (factor_base[i].p > BLOCK_SIZE)
            break;
    }

    r0 = i;

    for (b = 1; b <= qs_inf->sieve_size / BLOCK_SIZE; b++)
    {
        B = sieve + b * BLOCK_SIZE;

        for (pind = qs_inf->small_primes; pind < r0; pind++)
        {
            if (soln2[pind] == 0)
                continue;

            p = factor_base[pind].p;
            size = factor_base[pind].size;
            d1 = xr2[pind];
            d2 = p - d1;
            Bp = (B - xr2[pind]);

            for (x = sieve + xr1[pind]; x <= Bp; )
            {
                (*x) += size, x += d1;
                (*x) += size, x += d2;
            }

            if (x <= B)
            {
                (*x) += size, x += d1;
                xr2[pind] = d2;
            }
            else { xr2[pind] = d1; }

            xr1[pind] = (x - sieve);
        }
    }

    for (b = 1; b <= qs_inf->sieve_size / BLOCK_SIZE; b++)
    {
        B = sieve + b * BLOCK_SIZE;

        for (pind = r0; pind < num_primes; pind++)
        {
            p = factor_base[pind].p;

            if (soln2[pind] == 0)
                continue;

            size = factor_base[pind].size;

            if ((x = sieve + xr1[pind]) <= B)
            {
                (*x) += size;
                x += xr2[pind];

                if (x <= B)
                {
                    (*x) += size;
                    x += p - xr2[pind];
                }
                else { xr2[pind] = p - xr2[pind]; }

                xr1[pind] = (x - sieve);
            }
            else { xr1[pind] = (x - sieve); }
        }
    }
}

/* check position 'i' in sieve array for smoothness */

slong qsieve_evaluate_candidate(qs_t qs_inf, slong i, unsigned char * sieve)
{
   slong bits, exp, extra_bits;
   mp_limb_t modp, prime;
   slong num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   fac_t * factor = qs_inf->factor;
   mp_limb_t * soln1 = qs_inf->soln1;
   mp_limb_t * soln2 = qs_inf->soln2;
   mp_limb_t * A_ind = qs_inf->A_ind;
   slong * small = qs_inf->small;
   mp_limb_t pinv;
   slong num_factors = 0;
   slong relations = 0;
   slong j, k;

   fmpz_t X, Y, res, p;
   fmpz_init(X);
   fmpz_init(Y);
   fmpz_init(res);
   fmpz_init(p);

   fmpz_set_ui(X, i - qs_inf->sieve_size / 2); /* X */

   fmpz_mul(Y, X, qs_inf->A);
   fmpz_add(Y, Y, qs_inf->B); /* Y = AX+B */
   fmpz_add(res, Y, qs_inf->B);

   fmpz_mul(res, res, X);
   fmpz_add(res, res, qs_inf->C); /* res = AX^2 + 2BX + C */

   bits = FLINT_ABS(fmpz_bits(res));
   bits -= BITS_ADJUST;
   extra_bits = 0;

   if (factor_base[0].p != 1) /* divide out powers of the multiplier */
   {
      fmpz_set_ui(p, factor_base[0].p);
      exp = fmpz_remove(res, res, p);
      if (exp) extra_bits += exp*qs_inf->factor_base[0].size;
      small[0] = exp;
   }
   else small[0] = 0;

   fmpz_set_ui(p, 2); /* divide out by powers of 2 */
   exp = fmpz_remove(res, res, p);

   extra_bits += exp;
   small[1] = exp;

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
            if ((modp == soln1[j]) || (modp == soln2[j]))
            {
               fmpz_set_ui(p, prime);
               exp = fmpz_remove(res, res, p);
               if (exp)
               {
                  extra_bits += qs_inf->factor_base[j].size;
                  factor[num_factors].ind = j;
                  factor[num_factors++].exp = exp;
               }
            }
         }
         else
         {
            fmpz_set_ui(p, prime);
            exp = fmpz_remove(res, res, p);
            factor[num_factors].ind = j;
            factor[num_factors++].exp = exp + 1;
         }
      }

      if (fmpz_cmp_ui(res, 1) == 0 || fmpz_cmp_si(res, -1) == 0) /* We've found a relation */
      {
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

         qs_inf->num_factors = num_factors;

         qsieve_write_to_file(qs_inf, UWORD(1), Y);

         qs_inf->full_relation++;

        /* relations += qsieve_insert_relation(qs_inf, Y); */  /* Insert the relation in the matrix */

       /*  if (qs_inf->num_relations >= qs_inf->buffer_size)
         {
            flint_printf("Error: too many duplicate relations!\n");
            flint_printf("s = %wd, bits = %wd\n", qs_inf->s, qs_inf->bits);
            abort();
         } */
      }
      else
      {
          fmpz_abs(res, res);

          if (fmpz_bits(res) <= FLINT_BITS)
          {
              prime = fmpz_get_ui(res);

              if (prime > factor_base[qs_inf->q_idx].p && prime < 64 * factor_base[qs_inf->num_primes - 1].p && n_is_prime(prime))
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

                  qs_inf->num_factors = num_factors;

                  /* store this partial in file */

                  qsieve_write_to_file(qs_inf, prime, Y, i);

                  qs_inf->edges++;

                  qsieve_add_to_hashtable(qs_inf, prime);
              }
          }
      }

      goto cleanup;
   }

cleanup:
   fmpz_clear(X);
   fmpz_clear(Y);
   fmpz_clear(res);
   fmpz_clear(p);

   return relations;
}


/* scan sieve array for possible candidate for smooth relations */

slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve)
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
            {
                rels += qsieve_evaluate_candidate(qs_inf, i, sieve);
            }

            i++;
        }

        j++;
    }

    return rels;
}

/* procedure to call polynomial initialization and sieving procedure */

slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve)
{
    slong relation = 0;

    while (qs_inf->columns < qs_inf->num_primes + qs_inf->extra_rels)
    {
        qsieve_compute_C(qs_inf);
        qsieve_do_sieving(qs_inf, sieve);

        relation += qsieve_evaluate_sieve(qs_inf, sieve);

        if (qs_inf->curr_poly == (1 << qs_inf->s))
            break;

        qsieve_init_poly_next(qs_inf);
    }

    return relation;
}
