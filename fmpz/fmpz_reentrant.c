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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"

fmpz _fmpz_new_mpz(void)
{
   __mpz_struct * mpz_ptr = (__mpz_struct *) malloc(sizeof(__mpz_struct));
   mpz_init(mpz_ptr);
   return PTR_TO_COEFF(mpz_ptr);
}

void _fmpz_clear_mpz(fmpz f)
{
   mpz_clear(COEFF_TO_PTR(f));	
}

__mpz_struct * _fmpz_promote(fmpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) *f = _fmpz_new_mpz(); // f is small so promote it first
	// if f is large already, just return the pointer
      
   return COEFF_TO_PTR(*f);
}

__mpz_struct * _fmpz_promote_val(fmpz_t f)
{
   fmpz c = *f;
	if (!COEFF_IS_MPZ(c)) // f is small so promote it
	{
	   *f = _fmpz_new_mpz();
	   __mpz_struct * mpz_ptr = COEFF_TO_PTR(*f);
		mpz_set_si(mpz_ptr, c);
		return mpz_ptr;
	} else // f is large already, just return the pointer
      return COEFF_TO_PTR(*f);
}

void _fmpz_demote_val(fmpz_t f)
{
   __mpz_struct * mpz_ptr = COEFF_TO_PTR(*f);

	long size = mpz_ptr->_mp_size;
	
	if (size == 0L) // value is zero
	{
		_fmpz_clear_mpz(*f);
		*f = 0;
	} else if (size == 1L) // value is positive and 1 limb
	{
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= (ulong) COEFF_MAX) 
		{
			_fmpz_clear_mpz(*f);
			*f = (fmpz) uval;
		}
	} else if (size == -1L) // value is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= (ulong) COEFF_MAX) 
		{
			_fmpz_clear_mpz(*f);
			*f = (fmpz) -uval;
		}
	}
	// don't do anything if value has to be multi precision
}


