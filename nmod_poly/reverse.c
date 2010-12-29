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
#include "nmod_poly.h"

void nmod_poly_reverse(nmod_poly_t output, nmod_poly_t input, long m)
{
    long i, min;
    mp_limb_t temp;
      
    nmod_poly_fit_length(output, m);
    
    if (input != output)
    {
        min = FLINT_MIN(m, input->length);
        
        for (i = 0; i < min; i++)
            output->coeffs[m - i - 1] = input->coeffs[i];

        for ( ; i < m; i++)
            output->coeffs[m - i - 1] = 0L;
   } else
   {
      for (i = 0; i < m/2; i++)
      {
         temp = i < input->length ? input->coeffs[i] : 0;
         
         if (m - i - 1 < input->length)
            input->coeffs[i] = input->coeffs[m - i - 1];
         else
            input->coeffs[i] = 0;

         input->coeffs[m - i - 1] = temp;
      }
      if (m & 1 && i >= input->length) input->coeffs[i] = 0;
   }

   output->length = m;
   _nmod_poly_normalise(output);
}