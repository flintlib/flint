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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

mp_limb_t
nmod_poly_factor_with_berlekamp(nmod_poly_factor_t result, nmod_poly_t input)
{
    nmod_poly_t monic_input;
    nmod_poly_factor_t sqfree_factors, factors;
    mp_limb_t leading_coeff;
    long i, len;

    len = input->length;

    if (len < 2)
    {
        if (len == 0)
            return 0;
        else
            return input->coeffs[0];
    }

    leading_coeff = nmod_poly_get_coeff_ui(input, nmod_poly_degree(input));

    nmod_poly_init_preinv(monic_input, input->mod.n, input->mod.ninv);
    nmod_poly_make_monic(monic_input, input);

    if (len == 2)
    {
        nmod_poly_factor_insert(result, monic_input, 1);
        nmod_poly_clear(monic_input);
        return input->coeffs[1];
    }

    nmod_poly_factor_init(sqfree_factors);
    nmod_poly_factor_squarefree(sqfree_factors, monic_input);
    nmod_poly_clear(monic_input);

    /* Run Berlekamp on each of the square-free factors */
    for (i = 0; i < sqfree_factors->num_factors; i++)
    {
        nmod_poly_factor_init(factors);
        nmod_poly_factor_berlekamp(factors, sqfree_factors->factors[i]);
        nmod_poly_factor_pow(factors, sqfree_factors->exponents[i]);
        nmod_poly_factor_concat(result, factors);
        nmod_poly_factor_clear(factors);
    }

    nmod_poly_factor_clear(sqfree_factors);
    return leading_coeff;   
}
