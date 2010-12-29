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

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_poly.h"

void nmod_poly_make_monic(nmod_poly_t output, nmod_poly_t input)
{
    mp_limb_t inv;
    
    if (input->length == 0)
    {
        printf("Exception: division by zero in nmod_poly_invert\n");
        abort();
    }

    inv = n_invmod(input->coeffs[input->length - 1], input->mod.n);
    nmod_poly_scalar_mul(output, input, inv);
}