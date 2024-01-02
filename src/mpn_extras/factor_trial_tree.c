/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
#include <pthread.h>

static pthread_once_t _factor_trial_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t _factor_trial_lock;
#endif

FLINT_TLS_PREFIX mp_ptr _factor_trial_tree[16 - (FLINT_BITS/32)];
FLINT_TLS_PREFIX int _factor_trial_tree_initialised = 0;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
void _tree_mutex_init(void)
{
    pthread_mutex_init(&_factor_trial_lock, NULL);
}
#endif

void _cleanup_trial_tree(void)
{
    slong i;

    for (i = 0; i < 13 - (FLINT_BITS/32); i++)
	flint_free(_factor_trial_tree[i]);

    _factor_trial_tree_initialised = 0;
}

void
_factor_trial_tree_init(void)
{
    slong i, j, k, m, n;
    const mp_limb_t * primes;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&_factor_trial_initialised, _tree_mutex_init);
    pthread_mutex_lock(&_factor_trial_lock);
#endif

    if (!_factor_trial_tree_initialised)
    {
	    primes = n_primes_arr_readonly(3512);

	    flint_register_cleanup_function(_cleanup_trial_tree);

        /*
	        Initialise mpn's in tree, each of which is thought of as an array
	        of products of primes
	        Note there are 3512 primes less than 32768
	    */
        for (i = 0; i < 13 - (FLINT_BITS/32); i++)
	    {
	        _factor_trial_tree[i] = (mp_ptr)
		        flint_malloc(4096/(FLINT_BITS/16)*sizeof(mp_limb_t));
        }

	    /* initialise products in first layer of tree */
	    for (i = 0, j = 0; i < 3512; i+=(FLINT_BITS/16), j++)
        {
#if FLINT64
	        _factor_trial_tree[0][j] =
		        primes[i]*primes[i + 1]*primes[i + 2]*primes[i + 3];
#else
	        _factor_trial_tree[0][j] =
		        primes[i]*primes[i + 1];
#endif
	    }

	    m = 0; /* level in tree that has been computed already */
	    n = 1; /* number of words per entry in that level */
        k = 3512/(FLINT_BITS/16); /* number of entries on that level */

        /* compute remaining levels of tree */
	    for ( ; m < 12 - (FLINT_BITS/32); m++, n*=2, k=(k+1)/2)
        {
            /* multiply entries in pairs */
	        for (i = 0, j = 0; j < k/2; j++, i+=(2*n))
            {
                flint_mpn_mul_n(_factor_trial_tree[m + 1] + i,
		        _factor_trial_tree[m] + i,
		        _factor_trial_tree[m] + i + n, n);
	        }

	        if ((k % 2) == 1) /* copy across last entry if k is odd */
	        {
	            mpn_copyi(_factor_trial_tree[m + 1] + i,
	                                         _factor_trial_tree[m] + i, n);
                flint_mpn_zero(_factor_trial_tree[m + 1] + i + n, n);
	        }
	    }

        _factor_trial_tree_initialised = 1;
    }

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&_factor_trial_lock);
#endif
}

int flint_mpn_factor_trial_tree(slong * factors,
		                mp_srcptr x, mp_size_t xsize, slong num_primes)
{
    slong i, j, m, n, n2, nmax;
    const mp_limb_t * primes;
    mp_ptr gtemp; /* temporary space for recursive gcd's */
    mp_ptr temp; /* temporary space for flint_mpn_gcd_full2 */
    slong rlimbs[13 - (FLINT_BITS/32)]; /* gcd lengths at each level */
    slong idx[13 - (FLINT_BITS/32)]; /* indexes of entries at each level */
    slong offset; /* offset into gtemp */
    slong valid_level; /* last level with valid gcd computed */
    int numfacs = 0; /* number of prime factors found */
    int gcd1; /* whether we hit gcd = 1, can prune this branch of tree */

    _factor_trial_tree_init();

    primes = n_primes_arr_readonly(num_primes);

    gtemp = (mp_ptr)
	    flint_malloc((3*4096/(FLINT_BITS/16) + xsize)*sizeof(mp_limb_t));
    temp = gtemp + 2*4096/(FLINT_BITS/16);

    /* compute gcd of x with top level in tree */
    m  = FLINT_MAX(FLINT_BIT_COUNT(num_primes) - (FLINT_BITS/32), 0); /* top level in tree */
    n  = 4096/(FLINT_BITS/16); /* number of words of integer in tree */
    for (i = 12 - (FLINT_BITS/32); i > m; i--)
        n /= 2;
    n2 = n; /* number of words per entry at this level */
    nmax = n; /* save for later */

    MPN_NORM(_factor_trial_tree[m] + 0, n);

    if (n == 0) /* nothing to be done */
    {
        flint_free(gtemp);
        return 0;
    }

    rlimbs[m] = flint_mpn_gcd_full2(gtemp, x, xsize,
		                                   _factor_trial_tree[m] + 0, n, temp);

    if (rlimbs[m] == 1 && gtemp[0] == 1) /* gcd is 1, no factors found */
    {
        flint_free(gtemp);
	    return 0;
    }

    /* set indexes into arrays of entries at each level */
    for (j = 0; j < m; j++)
       idx[j] = -WORD(1);
    idx[m] = 0;

    valid_level = m; /* last level with valid gcd */

    /* compute gcd with each entry at bottom level of tree */
    for (i = 0; i < (num_primes + (FLINT_BITS/16) - 1)/(FLINT_BITS/16); i++)
    {
        gcd1 = 0;
	    n2 = nmax; /* number of words per entry */
        offset = 0; /* offset into gtemp */

	    /* compute gcd's down the tree for this i */
	    for (j = m; j >= 0; j--)
	    {
            /* must update indexes for current i */
	        if ((idx[j] & 1) != ((i >> j) & 1))
                idx[j]++;

            /* see if we need to compute new gcd for this level */
	        if (!gcd1 && (valid_level > j || (idx[j] & 1) != ((i >> j) & 1)))
            {
	            n = n2;
	            MPN_NORM(_factor_trial_tree[j] + idx[j]*n2, n);

                rlimbs[j] = flint_mpn_gcd_full2(gtemp + offset,
		             _factor_trial_tree[j] + idx[j]*n2, n,
		                               gtemp + offset - 2*n2, rlimbs[j + 1], temp);

                valid_level = j;

	            if (rlimbs[j] == 1 && gtemp[offset] == 1)
		            gcd1 = 1;
            }

	        offset += n2;
            n2/=2;
	    }

        if (!gcd1)
	    {
            for (j = 0; j < FLINT_BITS/16; j++)
            {

                /* check divisibility by primes with index 4*i + j */
	            if (flint_mpn_divisible_1_odd(x, xsize, primes[(FLINT_BITS/16)*i + j]))
	                factors[numfacs++] = (FLINT_BITS/16)*i + j;
	        }
	    }
    }

    flint_free(gtemp);

    return numfacs;
}

