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

    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"

void _nmod_poly_make_monic(mp_ptr output, 
                            mp_srcptr input, long len, nmod_t mod)
{
    mp_limb_t inv;
    
    inv = n_invmod(input[len - 1], mod.n);
    _nmod_vec_scalar_mul_nmod(output, input, len, inv, mod);
}

void nmod_poly_make_monic(nmod_poly_t output, const nmod_poly_t input)
{
    if (input->length == 0)
    {
        printf("Exception (nmod_poly_make_monic). Division by zero.\n");
        abort();
    }

    nmod_poly_fit_length(output, input->length);
    _nmod_poly_make_monic(output->coeffs, 
                            input->coeffs, input->length, input->mod);
    output->length = input->length;
}

