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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void nmod_poly_derivative(nmod_poly_t x_prime, nmod_poly_t x)
{
	long j;
    const long length = x->length;
	const mp_limb_t n = x->mod.n;
	const mp_limb_t ninv = x->mod.ninv;
	mp_limb_t k;
    
	if (length <= 1) 
	{
	   nmod_poly_zero(x_prime);
	   return;
    }

    nmod_poly_fit_length(x_prime, length - 1);	
	
	k = 1;
	for (j = 1; j < length; j++)
	{
		if (k <= 1) 
        {
            if (k == 0) x_prime->coeffs[j - 1] = 0L; 
		    else x_prime->coeffs[j - 1] = x->coeffs[j];
        } else 
            x_prime->coeffs[j - 1] = n_mulmod2_preinv(x->coeffs[j], k, n, ninv);
		
        if (++k == n) k = 0L;
	}

	x_prime->length = length - 1;
	_nmod_poly_normalise(x_prime); 
}