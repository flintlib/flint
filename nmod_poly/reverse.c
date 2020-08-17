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
#include "nmod_poly.h"

void _nmod_poly_reverse(mp_ptr output, mp_srcptr input, slong len, slong m)
{
    slong i, min;
    mp_limb_t temp;
      
    if (input != output)
    {
        min = FLINT_MIN(m, len);
        
        for (i = 0; i < min; i++)
            output[m - i - 1] = input[i];

        for ( ; i < m; i++)
            output[m - i - 1] = WORD(0);
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

void nmod_poly_reverse(nmod_poly_t output, const nmod_poly_t input, slong m)
{
    nmod_poly_fit_length(output, m);
     
    _nmod_poly_reverse(output->coeffs, input->coeffs, input->length, m);

    output->length = m;
    _nmod_poly_normalise(output);
}
