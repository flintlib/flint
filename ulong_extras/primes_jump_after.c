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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "ulong_extras.h"

void
n_primes_jump_after(n_primes_t iter, mp_limb_t n)
{
    if (n < iter->small_primes[iter->small_num - 1])
    {
        iter->small_i = n_prime_pi(n);
        iter->sieve_a = iter->sieve_b = iter->sieve_num = 0;
    }
    else
    {
        iter->small_i = iter->small_num;
        n_primes_sieve_range(iter, n + 1,
            n + 1 + FLINT_MIN(n, FLINT_SIEVE_SIZE) - 2);
    }
}
