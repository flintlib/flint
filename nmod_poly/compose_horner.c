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

void _nmod_poly_compose_horner(mp_ptr res, mp_srcptr poly1, 
                     long len1, mp_srcptr poly2, long len2, nmod_t mod)
{
    mp_limb_t * val1, * val2, t;
    long i, m, m2;
   
    if (len1 == 1) /* constant only */
    {
        res[0] = poly1[0];

        return;
    }

    if (len2 == 1) /* evaluate at constant */
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);

        return;
    }

    if (len1 == 2) /* linear poly not dealt with by general case */
    {
        t = poly1[0];
        _nmod_vec_scalar_mul_nmod(res, poly2, len2, poly1[1], mod);
        res[0] = n_addmod(res[0], t, mod.n);
       
        return;
    }

   	/* general case */
    m = len1 - 1;
	m2 = len2 - 1;

    val1 = nmod_vec_init(m*m2 + 1);
    val2 = nmod_vec_init(m*m2 + 1);

    /* initial c_m * poly2 + c_{m-1} */
    _nmod_vec_scalar_mul_nmod(val1, poly2, len2, poly1[m], mod);
    val1[0] = n_addmod(val1[0], poly1[m - 1], mod.n);
    
	m -= 2;
    i = 1;

	for ( ; m > 0; m--, i++) /* all but final val * poly2 + c_{m-1} */
	{
       _nmod_poly_mul(val2, val1, m2*i + 1, poly2, len2, mod);
       MP_PTR_SWAP(val1, val2);
       val1[0] = n_addmod(val1[0], poly1[m], mod.n);
	}

    /* final val * poly2 + c_0 */
    t = poly1[0];
    _nmod_poly_mul(res, val1, m2*(len1 - 2) + 1, poly2, len2, mod);
	res[0] = n_addmod(res[0], t, mod.n);

    nmod_vec_free(val1);
    nmod_vec_free(val2);
}

void nmod_poly_compose_horner(nmod_poly_t res, 
                       const nmod_poly_t poly1, const nmod_poly_t poly2)
{
	long len_out;
    
    if (poly1->length == 0) /* nothing to evaluate */ 
	{
	   nmod_poly_zero(res);

	   return;
	}

    if (poly2->length == 0) /* constant only */ 
	{
	   nmod_poly_set_coeff_ui(res, 0, poly1->coeffs[0]);
       nmod_poly_truncate(res, 1);

	   return;
	}

    len_out = (poly1->length - 1)*(poly2->length - 1) + 1;

	nmod_poly_fit_length(res, len_out);

    _nmod_poly_compose_horner(res->coeffs, poly1->coeffs, poly1->length, 
        poly2->coeffs, poly2->length, poly1->mod);

    res->length = len_out;

    _nmod_poly_normalise(res);

	return;
}
