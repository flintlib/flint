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

    Copyright (C) 2013 Marcin Bodych

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_smallest_primitive_root(mp_limb_t n)
{
    int i;
    mp_limb_t g;
    n_factor_t fac;
    mp_limb_signed_t exp;
    mp_limb_t phi;
	int ok;
	phi = n_euler_phi(n);
    n_factor_init(&fac);
    n_factor(&fac, phi, 1);
	g = 1;
	while( n > 2)
	{
		g++;
		if(g % 4 == 0 || g % 9 == 0 || g % 25 == 0) g++;
		ok = 1;
		if(n_powmod(g, phi, n) != 1) 
	    {
	 	    continue;
	    }
		for (i = fac.num-1; i >= 0; i--)
		{
			exp = phi / fac.p[i];
			if(n_powmod(g, exp, n) == 1) 
			{
				ok = 0;
				break;
			}
		}
		
		if(ok == 1)
			return g;
	}
    return g;
}
