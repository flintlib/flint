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

    Authored 2015 by Daniel S. Roche and A. Whitman Groves;
    US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"

void
fmpz_sparse_clear(fmpz_sparse_t poly)
{
  
  if(poly->coeffs)
  {
    slong i;
    for (i=0; i<poly->length; ++i) 
    {
      _fmpz_demote(poly->coeffs + i);
      _fmpz_demote(poly->expons + i);
    }
    flint_free(poly->coeffs);
    flint_free(poly->expons);
  }
}
