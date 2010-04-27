/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2010 William Hart
   
*****************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"

void mpfr_poly_realloc(mpfr_poly_t poly, const ulong alloc)
{
   ulong i;
   
   if (!alloc) // alloc == 0, clear up, reinitialise
   {
      mpfr_poly_clear(poly);
      mpfr_poly_init(poly, poly->prec);
	  return;
   }  
   
	if (poly->alloc) // realloc
	{
	   for (i = alloc; i < poly->alloc; i++)
	      mpfr_clear(poly->coeffs + i);

	   poly->coeffs = (mpfr *) realloc(poly->coeffs, alloc*sizeof(mpfr));
		
	   for (i = poly->alloc; i < alloc; i++)
			mpfr_init2(poly->coeffs + i, poly->prec);
	} else // nothing allocated already so do it now
	   poly->coeffs = _mpfr_vec_init(alloc, poly->prec);
   
   poly->alloc = alloc;  
}
