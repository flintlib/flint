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

    Copyright (C) 2015 Elena Sergeicheva
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_is_probabprime_miller_rabin(mp_limb_t n, mp_limb_t i)
{
	mp_limb_t g, j, rem, tmp, n_1, p, q, p_res;

	/* trivial */
	if (n == UWORD(2))
		return 1;
	if (n < UWORD(2) || (n & UWORD(1)) == 0 )
		return 0;

	/* check gcd(n, b) == 1 */
	if (i < 2) i = 2;
	for (; (g = n_gcd (n, i)) != 1; ++i)
		if (n > g) return 0;

	/* factor n-1 = q*2^p */
	n_1 = n-1;
	p_res = 0;
	tmp = n_1;
	while ((tmp & UWORD(1))==0)
	{
	  ++p_res;
	  tmp >>= 1; /* bisect */
	}
	p = p_res;
	q = tmp;

    if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
    {
    	/* calc b^q mod n; n is prime (or probabprime) in the case of 1 or n-1 */
	    rem = n_powmod(i, q, n);
		if (rem == 1 || rem == n_1)	return 1;

    	/* calc b^2q, b^4q, ... , b^((n-1)/2) mod n
    	   if any of them equals n-1, than n is prime (or probabprime) */
		for (j = 1; j < p; j++)
		{
			rem = n_powmod(rem, 2, n);
			if (rem == n_1)	return 1;
		}
    }
    else
    {

    	/* calc b^q mod n; n is prime (or probabprime) in the case of 1 or n-1 */
    	rem = n_powmod2(i, q, n);
    	if (rem == 1 || rem == n_1)	return 1;

    	/* calc b^2q, b^4q, ... , b^((n-1)/2) mod n
    	   if any of them equals n-1, than n is prime (or probabprime) */
    	for (j = 1; j < p; j++)
    	{
    		rem = n_powmod2(rem, 2, n);
    		if (rem == n_1) return 1;
    	}

    }

	return 0;
}
