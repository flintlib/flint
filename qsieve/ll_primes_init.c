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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"

prime_t * compute_factor_base(mp_limb_t * small_factor, qs_t qs_inf, long num_primes)
{
    mp_limb_t p, nmod, nmod2;
    mp_limb_t pinv;
    mp_limb_t k = qs_inf->k;
    long num = qs_inf->num_primes;
    long fb_prime = 2;
    prime_t * factor_base;
    int * sqrts;
    int kron;
    
    /* (re)allocate space for factor base */
    if (num == 0)
        factor_base = (prime_t *) flint_malloc(num_primes*sizeof(prime_t));
    else
        factor_base = (prime_t *) flint_realloc(qs_inf->factor_base, 
                                          num_primes*sizeof(prime_t));
    qs_inf->factor_base = factor_base;
    
    /* allocate space for square roots kn mod factor base primes */
    if (num == 0)
        sqrts = flint_malloc(sizeof(int)*num_primes);
    else
        sqrts = flint_realloc(qs_inf->sqrts, sizeof(int)*num_primes);
    qs_inf->sqrts = sqrts;

    qs_inf->num_primes = num_primes;

    /* compute the two limbs of kn */
    if (num == 0)
    {
        p = 2;
        num = 2;
    } else
        p = factor_base[num - 1].p;

    for (fb_prime = num; fb_prime < num_primes; ) /* leave space for k and 2 */
    {
        p = n_nextprime(p, 0);
        pinv = n_preinvert_limb(p);
        nmod = n_ll_mod_preinv(qs_inf->hi, qs_inf->lo, p, pinv); /* n mod p */
        if (nmod == 0) 
        {
            *small_factor = p;
            return factor_base;
        }
        
        nmod2 = n_mulmod2_preinv(nmod, k, p, pinv); /* kn mod p */
        if (nmod2 == 0) /* don't sieve with factors of multiplier */
            continue;
        
        nmod = nmod2; /* save nmod2 */

        kron = 1; /* n mod p is even, not handled by n_jacobi */
        while ((nmod2 % 2) == 0) 
        {
            if ((p % 8) == 3 || (p % 8) == 5) kron *= -1;
            nmod2 /= 2;
        }
        
        kron *= n_jacobi(nmod2, p); 
        if (kron == 1) /* kn is a quadratic residue mod p (and hence a FB prime) */
        {
            factor_base[fb_prime].p = p;
            factor_base[fb_prime].pinv = pinv;
            factor_base[fb_prime].size = FLINT_BIT_COUNT(p);
            sqrts[fb_prime] = n_sqrtmod(nmod, p);
            fb_prime++;
        }   
    }

    *small_factor = 0;
    return factor_base;
}

mp_limb_t qsieve_ll_primes_init(qs_t qs_inf)
{
    long num_primes;
    long i, s, min, fact, span;
    mp_limb_t fact_approx;
    fmpz_t temp;
    mp_limb_t k = qs_inf->k;
    mp_limb_t small_factor = 0;

    prime_t * factor_base;
    
    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < QS_LL_TUNE_SIZE; i++)
    {
        if (qsieve_ll_tune[i][0] > qs_inf->bits)
            break;
    }
    i--;
    
    qs_inf->sieve_bits = 32; /* number of bits to exceed in evaluate_sieve */
    qs_inf->sieve_size = qsieve_ll_tune[i][4]; /* size of sieve to use */
    qs_inf->small_primes = qsieve_ll_tune[i][3]; /* number of primes to not sieve with */
    num_primes = qsieve_ll_tune[i][2]; /* number of factor base primes */
    qs_inf->qsort_rels = qsieve_ll_tune[i][1]; /* number of relations to accumulate before sorting */
    
    qs_inf->num_primes = 0; /* start with 0 primes */
    factor_base = compute_factor_base(&small_factor, qs_inf, num_primes); /* build up FB */
    if (small_factor)
        return small_factor;

    /* figure out the number of factors of A and min, fact and span */
    s = qs_inf->bits/28 + 1; /* number of prime factors in A coeff */
   
    fmpz_init(temp);
    
    fmpz_mul_2exp(temp, qs_inf->kn, 1);
    fmpz_sqrt(temp, temp);
    fmpz_tdiv_q_ui(temp, temp, qs_inf->sieve_size);
    qs_inf->target_A = 2*fmpz_get_ui(temp);
   
    fmpz_root(temp, temp, s);
    fact_approx = fmpz_get_ui(temp);
   
    fmpz_clear(temp);

    fact = 2; 
    while (fact_approx >= factor_base[fact].p)
        fact++;
   
    while (1)
    {
        span = num_primes/s/s/2;
        if (span < 6*s) span = 6*s; /* make sure we have plenty of primes to choose from */

        min = fact - span/2;
        if (min < qs_inf->small_primes) 
            min = qs_inf->small_primes;

        fact = min + span/2;
        if (min + span <= num_primes - 2) /* we have enough primes */
            break;

        num_primes = (long) (1.1 * (double) num_primes);
        factor_base = compute_factor_base(&small_factor, qs_inf, num_primes); /* increase size of FB */
        if (small_factor)
            return small_factor;
    }
   
    qs_inf->s = s;
    qs_inf->min = min;
    qs_inf->fact = fact;
    qs_inf->span = span;

    /* these are used for the range of the middle factor of A when s is odd */
    qs_inf->mid = qs_inf->min + ((s - 1) * qs_inf->span)/(2 * s);
    qs_inf->high = qs_inf->mid + qs_inf->span / s;
   
#if (QS_DEBUG & 2)
    printf("Using %ld factor base primes\n", qs_inf->num_primes);
    printf("min = FB[%ld], span = %ld, number of A factors = %ld, target A = %ld\n", 
           min, span, s, qs_inf->target_A);
#endif
   
    /* consider k and 2 as factor base primes */
    factor_base[0].p = k;
    factor_base[0].pinv = n_preinvert_limb(k);
    factor_base[0].size = FLINT_BIT_COUNT(k);
    factor_base[1].p = 2;
    factor_base[0].size = 2;

    return 0;
}
