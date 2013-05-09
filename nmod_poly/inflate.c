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
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "nmod_poly.h"
#include "ulong_extras.h"

void
nmod_poly_inflate(nmod_poly_t result, const nmod_poly_t input, ulong inflation)
{
    if (input->length <= 1 || inflation == 1)
    {
        nmod_poly_set(result, input);
    }
    else if (inflation == 0)
    {
        mp_limb_t v = nmod_poly_evaluate_nmod(input, 1);
        nmod_poly_zero(result);
        nmod_poly_set_coeff_ui(result, 0, v);
    }
    else
    {
        len_t i, j, res_length = (input->length - 1) * inflation + 1;

        nmod_poly_fit_length(result, res_length);

        for (i = input->length - 1; i > 0; i--)
        {
            result->coeffs[i * inflation] = input->coeffs[i];
            for (j = i * inflation - 1; j > (i - 1) * inflation; j--)
                result->coeffs[j] = 0;
        }
        result->coeffs[0] = input->coeffs[0];
        result->length = res_length;
    }
}
