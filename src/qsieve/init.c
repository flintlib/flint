/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "qsieve.h"

#if (defined(__WIN32) && !defined(__CYGWIN__)) || defined(_MSC_VER)
#include <windows.h>
#endif

void qsieve_init_with_tune(qs_t qs_inf, const fmpz_t n, ulong ks_primes,
        slong fb_primes, slong small_primes, slong sieve_size, ulong sieve_bits)
{
    size_t fname_alloc_size;

#if (defined(__WIN32) && !defined(__CYGWIN__)) || defined(_MSC_VER)
    fname_alloc_size = MAX_PATH;
#else
    fname_alloc_size = 20;
#endif
    qs_inf->fname = (char *) flint_malloc(fname_alloc_size); /* space for filename */

    /* store n in struct */
    fmpz_init_set(qs_inf->n, n);

    /* determine the number of bits of n */
    qs_inf->bits =  fmpz_bits(n);

    /* store the tuning parameters in the struct */
    qs_inf->ks_primes    = ks_primes;    /* number of Knuth-Schroeppel primes */
    qs_inf->fb_primes    = fb_primes;    /* number of factor base primes */
    qs_inf->small_primes = small_primes; /* number of primes to not sieve with */
    qs_inf->sieve_size   = sieve_size;   /* size of sieve to use */

    /* split sieve_bits into threshold and fill (sieve values are biased so
       that the threshold is at least 64) */
    if (sieve_bits >= 64)
    {
        qs_inf->sieve_bits = sieve_bits;
        qs_inf->sieve_fill = 0;
    }
    else
    {
        qs_inf->sieve_bits = 64;
        qs_inf->sieve_fill = 64 - sieve_bits;
    }

    qs_inf->num_primes  = 0;
    qs_inf->num_relations = 0;
    qs_inf->full_relation = 0;
    qs_inf->num_cycles = 0;
    qs_inf->columns = 0;
    qs_inf->vertices = 0;
    qs_inf->components = 0;
    qs_inf->edges = 0;
#if QS_DEBUG
    qs_inf->poly_count = 0;
#endif

    fmpz_init(qs_inf->kn); /* initialise kn */

    qs_inf->factor_base = NULL;
    qs_inf->sqrts       = NULL;

    qs_inf->s = 0;
    qs_inf->low = 0;
    qs_inf->high = 0;
}

void qsieve_init(qs_t qs_inf, const fmpz_t n)
{
    slong i;
    flint_bitcnt_t bits = fmpz_bits(n);

    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < QS_TUNE_SIZE; i++)
    {
        if (qsieve_tune[i][0] > bits)
            break;
    }
    i--;

    /* initialise with the default tuning parameters for this bit-size */
    qsieve_init_with_tune(qs_inf, n,
        qsieve_tune[i][1],  /* ks_primes    */
        qsieve_tune[i][2],  /* fb_primes    */
        qsieve_tune[i][3],  /* small_primes */
        qsieve_tune[i][4],  /* sieve_size   */
        qsieve_tune[i][5]); /* sieve_bits   */
}
