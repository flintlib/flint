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
#include "nmod_poly.h"

void _nmod_poly_reverse(mp_ptr output, mp_srcptr input, long len, long m)
{
    long i, min;
    mp_limb_t temp;
      
    if (input != output)
    {
        min = FLINT_MIN(m, len);
        
        for (i = 0; i < min; i++)
            output[m - i - 1] = input[i];

        for ( ; i < m; i++)
            output[m - i - 1] = 0L;
    } else
    {
        for (i = 0; i < m/2; i++)
        {
            temp = i < len ? input[i] : 0;
         
            output[i] = m - i - 1 < len ? input[m - i - 1] : 0;

            output[m - i - 1] = temp;
        }
        if (m & 1 && i >= len) output[i] = 0;
    }
}

void nmod_poly_reverse(nmod_poly_t output, const nmod_poly_t input, long m)
{
    nmod_poly_fit_length(output, m);
     
    _nmod_poly_reverse(output->coeffs, input->coeffs, input->length, m);

    output->length = m;
    _nmod_poly_normalise(output);
}
