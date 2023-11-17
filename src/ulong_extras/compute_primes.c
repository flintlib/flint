/*
    Copyright (C) 2009 Tom Boothby
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
#include <pthread.h>

static pthread_once_t primes_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t primes_lock;
#endif

const unsigned int flint_primes_small[] =
{
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
    191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,
    281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
    389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
    491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,
    607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,
    719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
    829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,
    953,967,971,977,983,991,997,1009,1013,1019,1021
};


/* _flint_primes[i] holds an array of 2^i primes */
FLINT_TLS_PREFIX mp_limb_t * _flint_primes[FLINT_BITS];
FLINT_TLS_PREFIX double * _flint_prime_inverses[FLINT_BITS];
FLINT_TLS_PREFIX int _flint_primes_used = 0;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
void n_compute_primes_init()
{
   pthread_mutex_init(&primes_lock, NULL);
}
#endif

void
n_compute_primes(ulong num_primes)
{
    int i, m;
    ulong num_computed;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&primes_initialised, n_compute_primes_init);
    pthread_mutex_lock(&primes_lock);
#endif

    m = FLINT_CLOG2(num_primes);

    if (_flint_primes_used == 0)
        flint_register_cleanup_function(n_cleanup_primes);

    if (m >= _flint_primes_used)
    {
        n_primes_t iter;

        num_computed = UWORD(1) << m;
        _flint_primes[m] = flint_malloc(sizeof(mp_limb_t) * num_computed);
        _flint_prime_inverses[m] = flint_malloc(sizeof(double) * num_computed);

        n_primes_init(iter);
        for (i = 0; i < num_computed; i++)
        {
            _flint_primes[m][i] = n_primes_next(iter);
            _flint_prime_inverses[m][i] = n_precompute_inverse(_flint_primes[m][i]);
        }
        n_primes_clear(iter);

        /* copy to lower power-of-two slots */
        for (i = m - 1; i >= _flint_primes_used; i--)
        {
            _flint_primes[i] = _flint_primes[m];
            _flint_prime_inverses[i] = _flint_prime_inverses[m];
        }
        _flint_primes_used = m + 1;
    }

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&primes_lock);
#endif
}

