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
   Copyright (C) 2010 Daniel Woodhouse

*****************************************************************************/

/*
   This is testing get_period_sequence, a function developed for 
   Tom Coates and Alessio Corti at Imperial College
*/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mpoly.h"
#include "longlong.h"
#include "nmod_vec.h"
#include "fmpz.h"


int main(void)
{

#if FLINT_BITS == 64   
   fmpz_t zeroCoefficients[100], check;
   long coefficients[4] = {1,1,-1,1};
   ulong exponents[4]; 
   int i;

   ulong monomial = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 1 << 42);

   ulong primes3[20] = {
      1125899906842597,
      1125899906842589,
      1125899906842573,
      1125899906842553,
      1125899906842511,
      1125899906842507,
      1125899906842493,
      1125899906842463,
      1125899906842429,
      1125899906842391,
      1125899906842357,
      1125899906842283,
      1125899906842273,
      1125899906842247,
      1125899906842201,
      1125899906842177,
      1125899906842079,
      1125899906842033,
      1125899906842021,
      1125899906842013,
   };

   printf("get_period_sequence....");
   fflush(stdout);

   exponents[0] = (ulong) 0;
   exponents[1] = (ulong) 2 + ((ulong) 1 << 21) + ((ulong) 1 << 42);
   exponents[2] = (ulong) 1 + ((ulong) 2 << 21) + ((ulong) 1 << 42);
   exponents[3] = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 2 << 42);
   
   get_period_sequence(zeroCoefficients, coefficients, exponents, 4, monomial, 100, primes3, 10, 3);

   fmpz_init(check);
   fmpz_set_str(check, "-1612207508215775948685323966297082670959348818240567745024", 10);

   if (!(fmpz_equal(zeroCoefficients[99], check))){
      fmpz_print(zeroCoefficients[99]);
      printf(" FAIL\n");
   }

   for (i = 0; i < 100; i++)
      fmpz_clear(zeroCoefficients[i]);
#endif
 
   printf("PASS\n");
   return 0;
}
