/*
    Copyright (C) 2012 Thomas M. DuBuisson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fmpz-conversions.h"
#ifdef LONGSLONG
#define flint_mpz_cmp_si mpz_cmp_si
#else
#include "gmpcompat.h"
#endif

#ifdef FLINT_USES_POPCNT
static __inline__ flint_bitcnt_t shortCount(slong val)
{
#ifdef LONGSLONG
    return __builtin_popcountl(val);
#else
    return __builtin_popcountll(val);  
#endif
}
#else
/* A naive implementation if neither your processor nor your compiler want to
 * do the work. */
static __inline__ flint_bitcnt_t shortCount(slong val)
{
    flint_bitcnt_t cnt;
    for (cnt=0; val; val >>= 1)
        cnt += val & 1;
    return cnt;
}
#endif

flint_bitcnt_t fmpz_popcnt(const fmpz_t c)
{
    fmpz c1 = *c;

    if (c1 < 0)
        return 0;
    else if (!COEFF_IS_MPZ(c1))
        return shortCount(c1);
    else
    {
        mpz_mock_ptr t = COEFF_TO_PTR(c1);
        if (flint_mpz_cmp_si((mpz_ptr) t, 0) < 0)
            return 0;
        else
            return mpz_popcount((mpz_ptr) t);
    }
}
