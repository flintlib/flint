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

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int fmpz_cmp(const fmpz_t f, const fmpz_t g)
{
	if (f == g) return 0; // aliased inputs
	
	int sign;
   
   if (!COEFF_IS_MPZ(*f)) 
	{
		if (!COEFF_IS_MPZ(*g)) 
		{
         if (*f < *g) return -1;
			else return (*f > *g);
		} else // f is small, g is large
      {
         sign = mpz_sgn(COEFF_TO_PTR(*g));
         if ((*f >= 0) && (sign < 0)) return 1;
         return -sign;
      }
	} else 
   {
      if (!COEFF_IS_MPZ(*g)) // f is large, and g is small
      {
         sign = mpz_sgn(COEFF_TO_PTR(*f));
         if ((*g >= 0) && (sign < 0)) return -1;
         return sign; 
      } else 
         return mpz_cmp(COEFF_TO_PTR(*f), COEFF_TO_PTR(*g)); 
   }
}