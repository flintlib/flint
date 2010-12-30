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

/* helper: add constant to constant coefficient of p */
void __nmod_poly_add_ui(nmod_poly_t p, ulong c)
{
    if (p->length) /* add to existing coeff */
	    p->coeffs[0] = n_addmod(p->coeffs[0], c, p->mod.n);
	else 
        nmod_poly_set_coeff_ui(p, 0, c);

    _nmod_poly_normalise(p);
}

void nmod_poly_compose_horner(nmod_poly_t res, 
                       const nmod_poly_t poly1, const nmod_poly_t poly2)
{
	long m;
    mp_limb_t t;
    nmod_poly_t val;
	
    if (poly1->length == 0) /* nothing to evaluate */ 
	{
	   nmod_poly_zero(res);

	   return;
	}

	if (poly1->length == 1 || poly2->length == 0) /* constant only */
	{
	   nmod_poly_set_coeff_ui(res, 0, poly1->coeffs[0]);
       nmod_poly_truncate(res, 1);

	   return;
	}

    if (poly2->length == 1) /* evaluate at constant */
    {
        t = nmod_poly_evaluate(poly1, poly2->coeffs[0]);

        nmod_poly_set_coeff_ui(res, 0, t);
        nmod_poly_truncate(res, 1);

        return;
    }

	if (poly1->length == 2) /* linear poly not dealt with by general case */
	{
		t = poly1->coeffs[0];
		
        nmod_poly_scalar_mul(res, poly2, poly1->coeffs[1]);
        __nmod_poly_add_ui(res, t);

		return;
	}

	/* general case */
    m = poly1->length - 1;
	
    nmod_poly_init(val, poly1->mod.n);

	/* initial c_m * poly1 + c_{m-1} */
    nmod_poly_scalar_mul(val, poly2, poly1->coeffs[m]);
    __nmod_poly_add_ui(val, poly1->coeffs[m - 1]);
	
	m -= 2;

	for ( ; m > 0; m--) /* all but final val * poly1 + c_{m-1} */
	{
       nmod_poly_mul(val, val, poly2);
       __nmod_poly_add_ui(val, poly1->coeffs[m]);
	}

    /* final val * poly1 + c_0 */
    t = poly1->coeffs[0];
    nmod_poly_mul(res, val, poly2);
	__nmod_poly_add_ui(res, t);

    nmod_poly_clear(val);

	return;
}