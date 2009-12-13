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

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_randtest(fmpz_poly_t f, ulong length, mp_bitcnt_t bits_in)
{
   ulong i;
   fmpz_poly_fit_length(f, length);

   for (i = 0; i < length; i++)
      fmpz_randtest(f->coeffs + i, bits_in);
   
   _fmpz_poly_normalise(f);
}

void fmpz_poly_randtest_not_zero(fmpz_poly_t f, ulong length, mp_bitcnt_t bits_in)
{
   if ((bits_in == 0) || (length == 0))
   {
      printf("Exception : 0 passed to fmpz_poly_randtest_not_zero\n");
      abort();
   }

   fmpz_poly_randtest(f, length, bits_in);
   if (f->length == 0) 
   {
      fmpz_set_ui(f->coeffs, 1UL);
      f->length = 1;
   }
}