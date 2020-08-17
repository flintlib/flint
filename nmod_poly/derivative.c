/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void _nmod_poly_derivative(mp_ptr x_prime, mp_srcptr x, slong len, nmod_t mod)
{
	slong j;
    mp_limb_t k = 1;

	for (j = 1; j < len; j++)
	{
		if (k <= 1) 
            x_prime[j - 1] = k == 0 ? WORD(0) : x[j];     
        else 
            x_prime[j - 1] = n_mulmod2_preinv(x[j], k, mod.n, mod.ninv);
		
        if (++k == mod.n) k = WORD(0);
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

