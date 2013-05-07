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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

const unsigned char FLINT_PRIME_PI_ODD_LOOKUP[] =
{
  0,2,3,4,4,5,6,6,7,8,8,9,9,9,10,11,11,11,12,12,13,14,14,15,15,15,16,16,16,17,
  18,18,18,19,19,20,21,21,21,22,22,23,23,23,24,24,24,24,25,25,26,27,27,28,29,
  29,30,30,30,30,30,30,30,31,31,32,32,32,33,34,34,34,34,34,35,36,36,36,37,37,
  37,38,38,39,39,39,40,40,40,41,42,42,42,42,42,43,44,44,45,46,46,46,46,46,46,
  47,47,47,47,47,47,48,48,49,50,50,51,51,51,52,53,53,53,53,53,54,54,54,55,55,
  55,56,56,56,57,58,58,58,59,59,60,61,61,61,61,61,62,62,62,62,62,62,62,63,63
};


ulong n_prime_pi(mp_limb_t n)
{
    ulong low, mid, high;

    if (n < FLINT_PRIME_PI_ODD_LOOKUP_CUTOFF)
    {
        if (n < 3)
            return (n == 2);
        return FLINT_PRIME_PI_ODD_LOOKUP[(n-1)/2];
    }

    n_prime_pi_bounds(&low, &high, n);
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
