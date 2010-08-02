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
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_realloc(fmpz_mpoly_t poly, ulong alloc)
{
   if (alloc == 0) // clear up
   {
      fmpz_mpoly_clear(poly);
	  poly->alloc  = 0;
	  poly->length = 0;
      poly->coeffs = NULL;
      poly->exps   = NULL;
	  
	  return;
   }  
   
   if (poly->alloc != alloc)
   {
      if (poly->alloc) // realloc
	  {
         _fmpz_mpoly_truncate(poly, alloc);
         
		 poly->coeffs = (fmpz *) realloc(poly->coeffs, alloc*sizeof(fmpz));
		 poly->exps = (fmpz *) realloc(poly->exps, alloc*sizeof(fmpz));
		 if (alloc > poly->alloc)
		 {
		    mpn_zero((mp_ptr) (poly->coeffs + poly->alloc), alloc - poly->alloc);
            mpn_zero((mp_ptr) (poly->exps + poly->alloc), alloc - poly->alloc);
		 }
	  } else // nothing allocated already so do it now
	  {
         poly->coeffs = (fmpz *) calloc(alloc, sizeof(fmpz));
		 poly->exps = (fmpz *) calloc(alloc, sizeof(fmpz));
	  }
   }
   
   poly->alloc = alloc;
}
