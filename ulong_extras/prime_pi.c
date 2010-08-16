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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

#define PRIME_PI_LOOKUP_SIZE 64

char PRIME_PI_LOOKUP[PRIME_PI_LOOKUP_SIZE] =
{
    0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9,10,10,11,11,11,
    11,11,11,12,12,12,12,13,13,14,14,14,14,15,15,15,15,15,15,16,16,16,16,16,16,
    17,17,18,18,18
};

ulong n_prime_pi(ulong n)
{
    double logn;
    ulong low, mid, high;

    /* if (n < 2)
        return 0; */
    if (n < PRIME_PI_LOOKUP_SIZE)
        return PRIME_PI_LOOKUP[n];

    logn = n/log(n);
    high = 1.25506*logn + 1;
    low = logn - 1;

    n_compute_primes(high+1);

    while (low < high)
    {
        mid = (low + high) / 2;
        if (n < flint_primes[mid-1])
            high = mid;
        else
            low = mid + 1;
    }

    return low-1;
}
