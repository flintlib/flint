/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

/* mpn_random2 -- Generate random numbers with relatively long strings
   of ones and zeroes.  Suitable for border testing.

Copyright 1992-1994, 1996, 2000-2002, 2004, 2012 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */

#define mpn_incr_u(p,incr) \
  do { \
    mp_limb_t __x; \
    mp_ptr __p = (p); \
	__x = *__p + (incr); \
	*__p = __x; \
	if (__x < (incr)) \
	  while (++(*(++__p)) == 0);\
  } while (0)

void
flint_mpn_rrandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t nbits)
{
    flint_bitcnt_t bi;
    mp_limb_t ranm;		/* buffer for random bits */
    unsigned int cap_chunksize, chunksize;
    mp_size_t i;

    /* Set entire result to 111..1  */
    i = BITS_TO_LIMBS(nbits) - 1;
    rp[i] = UWORD_MAX >> (FLINT_BITS - (nbits % FLINT_BITS)) % FLINT_BITS;
    for (i = i - 1; i >= 0; i--)
        rp[i] = ~UWORD(0);

    ranm = n_randlimb(state);
    cap_chunksize = nbits / (ranm % 4 + 1);
    cap_chunksize += cap_chunksize == 0; /* make it at least 1 */

    bi = nbits;

    for (;;)
    {
        ranm = n_randlimb(state);
        chunksize = 1 + ranm % cap_chunksize;
        bi = (bi < chunksize) ? 0 : bi - chunksize;

        if (bi == 0)
            break;  /* low chunk is ...1 */

        rp[bi / FLINT_BITS] ^= UWORD(1) << (bi % FLINT_BITS);
        ranm = n_randlimb(state);
        chunksize = 1 + ranm % cap_chunksize;
        bi = (bi < chunksize) ? 0 : bi - chunksize;

        mpn_incr_u(rp + bi / FLINT_BITS, UWORD(1) << (bi % FLINT_BITS));

        if (bi == 0)
            break;      /* low chunk is ...0 */
    }
}

void flint_mpn_rrandom(mp_ptr rp, flint_rand_t state, mp_size_t n)
{
    mp_limb_t r = n_randlimb(state);

    if (r & (UWORD(1) << (FLINT_BITS - 1)))
        flint_mpn_rrandomb(rp, state, n * FLINT_BITS);
    else
        flint_mpn_rrandomb(rp, state, n * FLINT_BITS - r % FLINT_BITS);
}

void flint_mpn_urandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t n)
{
    slong i;
    slong nlimbs, nbits;

    nlimbs = (n + FLINT_BITS - 1) / FLINT_BITS;
    nbits = n % FLINT_BITS;

    for (i = 0; i < nlimbs; i++)
        rp[i] = _n_randlimb(state);

    if (nbits != 0)
        rp[nlimbs - 1] >>= (FLINT_BITS - nbits);
}

void flint_mpn_urandomm(mp_ptr rp, flint_rand_t state, mp_srcptr xp, mp_size_t xn)
{
    slong nbits;
    mp_limb_t lead;

    FLINT_ASSERT(xn > 0);
    FLINT_ASSERT(xp[xn - 1] != 0);

    if (xn == 1)
    {
        rp[0] = n_randint(state, xp[0]);
        return;
    }

    lead = xp[xn - 1];
    nbits = (xn - 1) * FLINT_BITS + FLINT_BIT_COUNT(lead);

    /* Special case for powers of two */
    if (flint_mpn_zero_p(xp, xn - 1) && (lead & (lead - 1)) == 0)
    {
        rp[xn - 1] = 0;  /* nbits - 1 means top limb might not be written */
        flint_mpn_urandomb(rp, state, nbits - 1);
    }
    else
    {
        /* Rejection sampling */
        do {
            flint_mpn_urandomb(rp, state, nbits);
        } while (mpn_cmp(rp, xp, xn) >= 0);
    }
}

