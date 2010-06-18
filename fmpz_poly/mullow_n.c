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

   Copyright (C) 2008, 2009 William Hart
   
*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void fmpz_poly_mullow_n(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong trunc)
{
   ulong limbs1 = fmpz_poly_max_limbs(poly1);
   ulong limbs2 = fmpz_poly_max_limbs(poly2);
   ulong max_limbs = FLINT_MAX(limbs1, limbs2);

   if (trunc <= 4)
   {
	   fmpz_poly_mullow_classical(res, poly1, poly2, trunc);
	   return;
   }
  
   if ((max_limbs > 3) && (trunc < 20))
      fmpz_poly_mullow_karatsuba_n(res, poly1, poly2, trunc);
   else
   {
	  fmpz_poly_mul_KS(res, poly1, poly2);
	  fmpz_poly_truncate(res, trunc);
   }
}
