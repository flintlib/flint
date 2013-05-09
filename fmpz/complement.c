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

   Copyright (C) 2012 Thomas M. DuBuisson

******************************************************************************/

#include "fmpz.h"

void fmpz_complement(fmpz_t r, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f)) /* f is small */
    {
	len_t res = ~(*f);
	fmpz_set_si(r, res);
    } else /* f is big */
    {
        if(r != f) { /* not aliased */
            __mpz_struct *ptr, *ptr2;
            ptr = _fmpz_promote(r);
            ptr2 = COEFF_TO_PTR(*f);
            mpz_com(ptr, ptr2);
            _fmpz_demote_val(r);
        } else { /* alaised */
            fmpz_t tmp;
            __mpz_struct *ptr, *ptr2;
            fmpz_init(tmp);
            ptr = _fmpz_promote(tmp);
            ptr2 = COEFF_TO_PTR(*f);
            mpz_com(ptr, ptr2);
            _fmpz_demote_val(tmp);
            fmpz_set(r,tmp);
            fmpz_clear(tmp);
        }
    }
}

