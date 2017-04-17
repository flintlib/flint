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

    Authored 2015 by A. Whitman Groves; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"
#include "fmpz_poly.h"

void fmpz_sparse_get_fmpz_poly(fmpz_poly_t out, const fmpz_sparse_t in)
{
    slong i, j;
    
    if(fmpz_sparse_is_zero(in))
    {
      fmpz_poly_zero(out);
      return;
    }
    else if(fmpz_sgn(in->expons) == -1)
    {
      fmpz_poly_zero(out);
      return;
    }
    else if(in->length == 1)
    {
      fmpz_poly_fit_length(out, 1);
      fmpz_set(out->coeffs, in->coeffs);
      out->length = 1;
      return;
    }
    else
    {
      i = 0;
      j = fmpz_get_si(in->expons);
      
      fmpz_poly_fit_length(out, j+1); 
      
      while (i < in->length && fmpz_cmp_si(in->expons + i, 0) >= 0)
      {
        fmpz_set(out->coeffs + fmpz_get_si(in->expons + i), in->coeffs + i);
        i++;
      }

      _fmpz_poly_set_length(out, j+1);
    }
}
