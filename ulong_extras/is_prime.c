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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int n_is_prime(mp_limb_t n)
{
#if !FLINT64
    return n_is_probabprime(n);
#else
    int isprime;
    if (n < 10000000000000000UL)
        return n_is_probabprime(n);

    if (!n_is_probabprime(n)) return 0;

    isprime = n_is_prime_pocklington(n, 100);
    if (isprime != -1) return isprime;

    return n_is_prime_pseudosquare(n);
#endif
}
