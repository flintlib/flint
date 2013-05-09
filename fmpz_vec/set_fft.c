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

    Copyright (C) 2008-2011 William Hart
    
******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fft.h"

void _fmpz_vec_set_fft(fmpz * coeffs_m, len_t length,
                          mp_limb_t ** coeffs_f, len_t limbs, len_t sign)
{
   len_t i, size;
   mp_limb_t * data;
   __mpz_struct * mpz_ptr;

   if (sign)
   {   
		for (i = 0; i < length; i++)
      {
         mpz_ptr = _fmpz_promote(coeffs_m);
         if (mpz_ptr->_mp_alloc < limbs) _mpz_realloc(mpz_ptr, limbs);
			data = mpz_ptr->_mp_d;
			
			if ((coeffs_f[i][limbs - 1] >> (FLINT_BITS - 1)) || coeffs_f[i][limbs])
         {
            mpn_neg_n(data, coeffs_f[i], limbs);
            mpn_add_1(data, data, limbs, 1L); 
            size = limbs;
			   while ((size) && (data[size - 1] == 0)) size--; /* normalise */
			   mpz_ptr->_mp_size = -size;
			   if (size >= -1L) _fmpz_demote_val(coeffs_m); /* coefficient may be small*/
         } else
         {
            flint_mpn_copyi(data, coeffs_f[i], limbs); 
			   size = limbs;
			   while ((size) && (data[size - 1] == 0L)) size--; /* normalise */
			   mpz_ptr->_mp_size = size;
			   if (size <= 1) _fmpz_demote_val(coeffs_m); /* coefficient may be small */
         }
         
			coeffs_m++;
      }
   } else 
   {
		for (i = 0; i < length; i++)
      {
         mpz_ptr = _fmpz_promote(coeffs_m);
         if (mpz_ptr->_mp_alloc < limbs) _mpz_realloc(mpz_ptr, limbs);
			data = mpz_ptr->_mp_d;
			flint_mpn_copyi(data, coeffs_f[i], limbs); 
			size = limbs;
			while ((size) && (data[size - 1] == 0L)) size--; /* normalise */
			mpz_ptr->_mp_size = size;
			if (size <= 1) _fmpz_demote_val(coeffs_m); /* coefficient may be small */

			coeffs_m++;
      }	
   }
}
