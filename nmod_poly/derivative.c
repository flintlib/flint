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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void _nmod_poly_derivative(mp_ptr x_prime, mp_srcptr x, len_t len, nmod_t mod)
{
	len_t j;
    mp_limb_t k = 1;

	for (j = 1; j < len; j++)
	{
		if (k <= 1) 
            x_prime[j - 1] = k == 0 ? 0L : x[j];     
        else 
            x_prime[j - 1] = n_mulmod2_preinv(x[j], k, mod.n, mod.ninv);
		
        if (++k == mod.n) k = 0L;
	}

}

void nmod_poly_derivative(nmod_poly_t x_prime, const nmod_poly_t x)
{
	if (x->length <= 1) 
	{
	   nmod_poly_zero(x_prime);
	   return;
    }

    nmod_poly_fit_length(x_prime, x->length - 1);	
	
    _nmod_poly_derivative(x_prime->coeffs, x->coeffs, x->length, x->mod);

	x_prime->length = x->length - 1;
	_nmod_poly_normalise(x_prime); 
}

