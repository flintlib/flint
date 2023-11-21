/*
    Copyright (C) 2006, 2011, 2016, 2020 William Hart
    Copyright (C) 2015 Nitin Kumar
    Copyright (C) 2020 Dan Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

#ifdef __GNUC__
# define memset __builtin_memset
#else
# include <string.h>
#endif

/*
    The actual sieving part of the quadratic sieve. This version is only run
    if the sieve block is small, i.e. less than 2*BLOCK_SIZE.
*/
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

   /* set entries in sieve to initial value and put a sentinel at the end */
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

      /* mark off entries in sieve, unrolled by 2 */
      while (bound - pos1 > 0)
      {
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
      }

      /* deal with final entries at the end, missed by unrolled loop */
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
   Second sieving routine. This version breaks sieve array into blocks then
   sieves each block, with the whole factor base at once.
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

    /* fill sieve with initial values and put sentinel at end */
    memset(sieve, qs_inf->sieve_fill, qs_inf->sieve_size + sizeof(ulong));
    sieve[qs_inf->sieve_size] = (char) 255;

    /*
       initial values for positions (which must be saved at the end of each
       sieve block in preparation for start of next sieve block)
    */
    for (i = 0; i < num_primes; i++)
    {
        posn1[i] = soln1[i];
        posn2[i] = soln2[i] - posn1[i];
    }

    /* break sieve into blocks of size BLOCK_SIZE */
    for (b = 1; b <= qs_inf->sieve_size / BLOCK_SIZE; b++)
    {
        /* end of current sieve block */
        B = sieve + b * BLOCK_SIZE;

        /*
            deal with small to medium sized primes first
            these hit sieve block multiple times, making unrolling worthwhile
        */
        for (pind = qs_inf->small_primes; pind < qs_inf->second_prime; pind++)
        {
            if (soln2[pind] == 0) /* skip primes dividing A */
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
            else
            {
                posn2[pind] = d1;
            }

            posn1[pind] = (pos - sieve);
        }

        /*
            now deal with larger primes which are likely to only hit sieve
            block once, if at all
        */
        for (pind = qs_inf->second_prime; pind < num_primes; pind++)
        {
            p = factor_base[pind].p;

            if (soln2[pind] == 0)
                continue;

            size = factor_base[pind].size;
            pos = sieve + posn1[pind];

            if (pos < B) /* hits the sieve block */
            {
                (*pos) += size;
                pos += posn2[pind];

                if (pos < B)
                {
                    (*pos) += size;
                    pos += p - posn2[pind];
                } else
                {
                    posn2[pind] = p - posn2[pind];
                }

                posn1[pind] = (pos - sieve);
            } else /* doesn't hit this sieve block, save posn for next block */
            {
                posn1[pind] = (pos - sieve);
            }
        }
    }
}

/*
    check position i in sieve array for smoothness
*/
slong qsieve_evaluate_candidate(qs_t qs_inf, ulong i, unsigned char * sieve, qs_poly_t poly)
{
   slong bits, exp, extra_bits;
   mp_limb_t modp, prime;
   slong num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   slong * small = poly->small; /* exponents of small primes and mult. */
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
   flint_printf("res = "); fmpz_print(res); flint_printf("\n");
   flint_printf("Poly: "); fmpz_print(qs_inf->A); flint_printf("*x^2 + 2*");
   fmpz_print(poly->B); flint_printf("*x + "); fmpz_print(C); flint_printf("\n");
   flint_printf("x = %wd\n", i - qs_inf->sieve_size / 2);
#endif

   sieve[i] -= qs_inf->sieve_fill; /* adjust sieve entry to number of bits */
   bits = FLINT_ABS(fmpz_bits(res)); /* compute bits of poly value */
   bits -= BITS_ADJUST; /* adjust for log approximations */
   extra_bits = 0; /* bits for mult. and small primes we didn't sieve with */

   if (factor_base[0].p != 1) /* divide out powers of the multiplier */
   {
      fmpz_set_ui(p, factor_base[0].p);
      exp = fmpz_remove(res, res, p);
      if (exp)
          extra_bits += exp*qs_inf->factor_base[0].size;
      small[0] = exp;
#if QS_DEBUG & 128
      if (exp)
          flint_printf("%ld^%ld ", factor_base[0].p, exp);
#endif
   }
   else
   {
       small[0] = 0;
   }

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
         if (exp)
             extra_bits += qs_inf->factor_base[j].size;
         small[j] = exp;
#if QS_DEBUG & 128
         if (exp)
            flint_printf("%ld^%ld ", prime, exp);
#endif
      }
      else
      {
         small[j] = 0;
      }
   }

   /* if we have a chance of being a candidate */
   if (extra_bits + sieve[i] > bits)
   {
      sieve[i] += extra_bits; /* add bits from small primes and mult. */

      /*
          pull out remaining primes
          we stop once we have reached the same number of bits as indicated in
          the sieve entry
      */
      for (j = qs_inf->small_primes; j < num_primes && extra_bits < sieve[i]; j++)
      {
         prime = factor_base[j].p;
         pinv = factor_base[j].pinv;
         modp = n_mod2_preinv(i, prime, pinv);

         if (soln2[j] != 0) /* not a prime dividing A */
         {
            if (modp == soln1[j] || modp == soln2[j])
            {
               fmpz_set_ui(p, prime);
               exp = fmpz_remove(res, res, p);
               extra_bits += qs_inf->factor_base[j].size;
               factor[num_factors].ind = j;
               factor[num_factors++].exp = exp;
#if QS_DEBUG & 128
                  flint_printf("%ld^%ld ", prime, exp);
#endif

            }
         } else /* prime dividing A */
         {
            fmpz_set_ui(p, prime);
            exp = fmpz_remove(res, res, p);
            factor[num_factors].ind = j;
            factor[num_factors++].exp = exp + 1; /* really factoring A*f(i) */
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
            flint_printf("%ld relations\n", qs_inf->full_relation);
#endif
         /* set sign amongst small factors */
         if (fmpz_cmp_si(res, -1) == 0)
            small[2] = 1;
         else
            small[2] = 0;

         for (k = 0; k < qs_inf->s; k++) /* Commit any outstanding A factors */
         {
            if (A_ind[k] >= j) /* check it is beyond where loop above ended */
            {
               factor[num_factors].ind = A_ind[k];
               factor[num_factors++].exp = 1;
            }
         }

         poly->num_factors = num_factors;

#if FLINT_USES_PTHREAD
         pthread_mutex_lock(&qs_inf->mutex);
#endif
	 qsieve_write_to_file(qs_inf, 1, Y, poly);

         qs_inf->full_relation++;

#if FLINT_USES_PTHREAD
         pthread_mutex_unlock(&qs_inf->mutex);
#endif
         relations++;
      } else /* not a relation, perhaps a partial? */
      {
          /* set sign */
          if (fmpz_sgn(res) < 0)
          {
              fmpz_neg(res, res);
              small[2] = 1;
          } else
              small[2] = 0;

          /* if we have a small cofactor (at most 30 bits) */
          if (fmpz_bits(res) <= 30)
          {
              prime = fmpz_get_ui(res);

              /*
                 a large prime is taken heuristically to be < 60 times largest
                 FB prime; skip values not coprime with multiplier, as this
                 will lead to factors of kn, not n
              */
              if (prime < 60*factor_base[qs_inf->num_primes - 1].p && n_gcd(prime, qs_inf->k) == 1)
              {
                  for (k = 0; k < qs_inf->s; k++)  /* commit any outstanding A factors */
                  {
                      if (A_ind[k] >= j) /* check beyond where loop above ended */
                      {
                          factor[num_factors].ind = A_ind[k];
                          factor[num_factors++].exp = 1;
                      }
                  }

                  poly->num_factors = num_factors;

#if FLINT_USES_PTHREAD
                  pthread_mutex_lock(&qs_inf->mutex);
#endif
                  /* store this partial in file */

                  qsieve_write_to_file(qs_inf, prime, Y, poly);

                  qs_inf->edges++;

                  qsieve_add_to_hashtable(qs_inf, prime);

#if FLINT_USES_PTHREAD
                  pthread_mutex_unlock(&qs_inf->mutex);
#endif

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

/*
    scan sieve array for possible candidates for smooth relations and process
    those relations to check for smoothness
*/
slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
    slong i = 0, j = 0;
    ulong * sieve2 = (ulong *) sieve;
    unsigned char bits = qs_inf->sieve_bits;
    slong rels = 0;

    while (j < qs_inf->sieve_size / sizeof(ulong))
    {
        /* scan 4 or 8 bytes at once for sieve entries over threshold */
#if FLINT64
        while ((sieve2[j] & UWORD(0xC0C0C0C0C0C0C0C0)) == 0)
#else
        while ((sieve2[j] & UWORD(0xC0C0C0C0)) == 0)
#endif
        {
            j++; /* advance to next word */
        }

        i = j * sizeof(ulong);

        /* check bytes individually in word */
        while (i < (j + 1) * sizeof(ulong) && i < qs_inf->sieve_size)
        {
            /* if we are over the threshold, check candidate for smoothness */
            if (sieve[i] > bits)
               rels += qsieve_evaluate_candidate(qs_inf, i, sieve, poly);

            i++;
        }

        j++;
    }

    return rels;
}

/* procedure to call polynomial initialization and sieving procedure */

typedef struct
{
    qs_s * inf;
    unsigned char * sieve;
    slong thread_idx;
    qs_poly_s * thread_poly;
    unsigned char * thread_sieve;
    slong rels;
}
_worker_arg_struct;

static void qsieve_collect_relations_worker(void * varg)
{
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    qs_s * qs_inf = arg->inf;
    qs_poly_s * thread_poly = arg->thread_poly;
    unsigned char * thread_sieve = arg->thread_sieve;
    slong j, iterations = (1 << (qs_inf->s - 1));

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&qs_inf->mutex);
#endif
        j = qs_inf->index_j;
        qs_inf->index_j = j + 1;
        if (j < iterations)
        {
           /* copy poly data for thread we are in */
            if (j > 0)
                qsieve_init_poly_next(qs_inf, j);
            qsieve_poly_copy(thread_poly, qs_inf);
        }
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&qs_inf->mutex);
#endif

        if (j >= iterations)
            return;

        if (qs_inf->sieve_size < 2*BLOCK_SIZE)
           qsieve_do_sieving(qs_inf, thread_sieve, thread_poly);
        else
           qsieve_do_sieving2(qs_inf, thread_sieve, thread_poly);

        arg->rels += qsieve_evaluate_sieve(qs_inf, thread_sieve, thread_poly);
    }
}


slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve)
{
    slong i;
    _worker_arg_struct * args;
    slong relations;
    const thread_pool_handle* handles = qs_inf->handles;
    slong num_handles = qs_inf->num_handles;

    args = (_worker_arg_struct *) flint_malloc((1 + num_handles)
                                                  *sizeof(_worker_arg_struct));

    qs_inf->index_j = 0;
    qsieve_init_poly_first(qs_inf);

    for (i = 0; i <= num_handles; i++)
    {
        args[i].inf = qs_inf;
        args[i].thread_idx = i;
        args[i].thread_poly = qs_inf->poly + i;
        args[i].thread_sieve = sieve + (qs_inf->sieve_size + sizeof(ulong) + 64)*i;
        args[i].rels = 0;
    }

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i], 0,
                                    qsieve_collect_relations_worker, &args[i]);
    }

    qsieve_collect_relations_worker(&args[num_handles]);

    relations = args[num_handles].rels;

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
        relations += args[i].rels;
    }

    flint_free(args);

    return relations;
}
