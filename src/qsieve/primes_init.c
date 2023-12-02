/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

prime_t *
compute_factor_base(mp_limb_t * small_factor, qs_t qs_inf, slong num_primes)
{
    mp_limb_t p, nmod, nmod2;
    mp_limb_t pinv;
    mp_limb_t k = qs_inf->k;
    slong num = qs_inf->num_primes;
    slong fb_prime = 2;
    prime_t * factor_base = NULL;
    int * sqrts;
    int kron;
    n_primes_t iter;

    /* (re)allocate space for factor base */

    factor_base = (prime_t *) flint_realloc(qs_inf->factor_base,
                                          num_primes*sizeof(prime_t));
    qs_inf->factor_base = factor_base;

    /* allocate space for square roots kn mod factor base primes */
    sqrts = flint_realloc(qs_inf->sqrts, sizeof(int)*num_primes);

    qs_inf->sqrts = sqrts;

    /* compute the last prime in factor base */
    p = num == 0 ? 2 : factor_base[num - 1].p;
    if (num == 0)
       num = 3;

    n_primes_init(iter);
    n_primes_jump_after(iter, p);

    /* factor base already contains num primes */
    for (fb_prime = num; fb_prime < num_primes; )
    {
        p = n_primes_next(iter);
        pinv = n_preinvert_limb(p);
        nmod = fmpz_fdiv_ui(qs_inf->n, p); /* n mod p */

        if (nmod == 0)
        {
            n_primes_clear(iter);
            *small_factor = p;
            return factor_base;
        }

        nmod2 = n_mulmod2_preinv(nmod, k, p, pinv); /* kn mod p */

        if (nmod2 == 0) /* don't sieve with factors of multiplier */
            continue;

        nmod = nmod2; /* save nmod2 */
        kron = 1; /* n mod p is even, not handled by n_jacobi */

        while (nmod2 % 2 == 0)
        {
            if (p % 8 == 3 || p % 8 == 5) kron *= -1;
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

    n_primes_clear(iter);

    *small_factor = 0;

    return factor_base;
}

mp_limb_t qsieve_primes_init(qs_t qs_inf)
{
    slong num_primes;
    slong i;
    mp_limb_t k = qs_inf->k;
    mp_limb_t small_factor = 0;
    slong bits;

    prime_t * factor_base;

    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < QS_TUNE_SIZE; i++)
    {
        if (qsieve_tune[i][0] > qs_inf->bits)
            break;
    }
    i--;

    num_primes = qsieve_tune[i][2]; /* number of factor base primes */

    if (num_primes < 3)
    {
       flint_throw(FLINT_ERROR, "Too few factor base primes\n");
    }

    qs_inf->sieve_size = qsieve_tune[i][4]; /* size of sieve to use */
    qs_inf->small_primes = qsieve_tune[i][3]; /* number of primes to not sieve with */

    bits = qsieve_tune[i][5];
    if (bits >= 64)
    {
       qs_inf->sieve_bits = bits;
       qs_inf->sieve_fill = 0;
    } else
    {
       qs_inf->sieve_bits = 64;
       qs_inf->sieve_fill = 64 - bits;
    }

    if (qs_inf->small_primes > num_primes)
    {
       flint_throw(FLINT_ERROR, "Too few factor base primes\n");
    }

    /* build up FB */
    factor_base = compute_factor_base(&small_factor, qs_inf, num_primes + qs_inf->ks_primes);
    if (small_factor)
        return small_factor;

    qs_inf->num_primes = num_primes;

    /* calculating optimal A coefficient size for hypercube */
    fmpz_init(qs_inf->target_A);
    fmpz_mul_2exp(qs_inf->target_A, qs_inf->kn, 1);
    fmpz_sqrt(qs_inf->target_A, qs_inf->target_A);
    fmpz_tdiv_q_ui(qs_inf->target_A, qs_inf->target_A, qs_inf->sieve_size/2);

    /* consider k and 2 and -1 as factor base primes */
    factor_base[0].p = k;
    factor_base[0].pinv = n_preinvert_limb(k);
    factor_base[0].size = FLINT_BIT_COUNT(k);
    factor_base[1].p = 2;
    factor_base[1].size = 2;
    factor_base[2].p = -1;

    return 0;
}

/* function to call for incrementing number of factor base prime by 'delta' */

mp_limb_t qsieve_primes_increment(qs_t qs_inf, mp_limb_t delta)
{
    slong num_primes = qs_inf->num_primes + delta;
    mp_limb_t small_factor = 0;

    compute_factor_base(&small_factor, qs_inf, num_primes + qs_inf->ks_primes);

    /* calculating optimal A coefficient size for hypercube */
    fmpz_init(qs_inf->target_A);
    fmpz_mul_2exp(qs_inf->target_A, qs_inf->kn, 1);
    fmpz_sqrt(qs_inf->target_A, qs_inf->target_A);
    fmpz_tdiv_q_ui(qs_inf->target_A, qs_inf->target_A, qs_inf->sieve_size/2);

    qs_inf->num_primes = num_primes;

    if (small_factor)
        return small_factor;

    return 0;
}
