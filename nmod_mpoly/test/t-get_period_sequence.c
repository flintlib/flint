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

// This is testing get_period_sequence, a function developed for Tom Coates and
// Alessio Corti at Imperial


#include <mpir.h>
#include <stdlib.h>
//#include <time.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mpoly.h"
#include "longlong.h"
#include "nmod_vec.h"
#include "fmpz.h"


int main(void)
{
   int result;
   printf("get_period_sequence....");
   fflush(stdout);

   fmpz_t zeroCoefficients[100];
   long coefficients[4] = {1,1,-1,1};
   ulong exponents[4]; 

   exponents[0] = (ulong) 0;
   exponents[1] = (ulong) 2 + ((ulong) 1 << 21) + ((ulong) 1 << 42);
   exponents[2] = (ulong) 1 + ((ulong) 2 << 21) + ((ulong) 1 << 42);
   exponents[3] = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 2 << 42);

   ulong monomial = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 1 << 42);

   ulong primes[30] = {104549, 104551, 104561, 104579, 104593, 104597, 104623, 104639, 104651, 104659, 
   104677, 104681, 104683, 104693, 104701, 104707, 104711, 104717, 104723, 104729,
   104417, 104459, 104471, 104473, 104479, 104491, 104513, 104527, 104537, 104543};

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

   /*
     clock_t start, end;
   double elapsed;
   start = clock();
   */

   get_period_sequence(zeroCoefficients, coefficients, exponents, 4, monomial, 100, primes3, 10, 3);
 
   /*
   end = clock();
   elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
   */
   

   fmpz_t check;
   fmpz_init(check);
   fmpz_set_str(check, "-1612207508215775948685323966297082670959348818240567745024", 10);

   if(!(fmpz_equal(zeroCoefficients[99], check))){
      fmpz_print(zeroCoefficients[99]);
      printf(" FAIL\n");
   }

   int i;
   /*
   for(i=0; i<100; i++){
      fmpz_print(zeroCoefficients[i]);
      printf("\n");
   }*/

   for(i=0; i<100; i++){
      fmpz_clear(zeroCoefficients[i]);
   }

   //check again with version2 of the function
   get_period_sequence2(zeroCoefficients, coefficients, exponents, 4, monomial, 100, primes3, 10, 3);
 
   /*
   end = clock();
   elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
   */
   

  
   
   fmpz_set_str(check, "-1612207508215775948685323966297082670959348818240567745024", 10);

   if(!(fmpz_equal(zeroCoefficients[99], check))){
      fmpz_print(zeroCoefficients[99]);
      printf(" FAIL\n");
   }

   fmpz_clear(check);
   /*
   for(i=0; i<100; i++){
      fmpz_print(zeroCoefficients[i]);
      printf("\n");
   }*/

   for(i=0; i<100; i++){
      fmpz_clear(zeroCoefficients[i]);
   }


   fmpz_t zeroCoefficients2[10];
   long coefficients2[9] = {1,1,1,1,2,1,1,1,1};
   ulong exponents2[9]; 

   exponents2[0] = (ulong) 0 + ((ulong) 2 << 21) + ((ulong) 0 << 42);
   exponents2[1] = (ulong) 0 + ((ulong) 1 << 21) + ((ulong) 1 << 42);
   exponents2[2] = (ulong) 1 + ((ulong) 2 << 21) + ((ulong) 1 << 42);
   exponents2[3] = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 2 << 42);
   exponents2[4] = (ulong) 1 + ((ulong) 2 << 21) + ((ulong) 2 << 42);
   exponents2[5] = (ulong) 1 + ((ulong) 0 << 21) + ((ulong) 3 << 42);
   exponents2[6] = (ulong) 2 + ((ulong) 1 << 21) + ((ulong) 3 << 42);
   exponents2[7] = (ulong) 1 + ((ulong) 2 << 21) + ((ulong) 3 << 42);
   exponents2[8] = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 4 << 42);

   ulong monomial2 = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 3 << 42);
 
   get_period_sequence(zeroCoefficients2, coefficients2, exponents2, 9, monomial2, 10, primes3, 10, 3);


   fmpz_set_str(check, "1348704", 10);

   if(!(fmpz_equal(zeroCoefficients2[9], check))){
      fmpz_print(zeroCoefficients2[9]);
      printf("\nFAIL\n");
   } 

   for(i=0; i<10; i++){
      fmpz_clear(zeroCoefficients2[i]);
   }
    
 /*
   fmpz_t zeroCoefficients3[100];
   long coefficients3[14] = {1,2,2,2,2,2,1,2,1,2,1,2,1,1};
   ulong exponents3[14]; 

   exponents3[0] = (ulong) 2 + ((ulong) 2 << 21) + ((ulong) 0 << 42);
   exponents3[1] = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 1 << 42);
   exponents3[2] = (ulong) 2 + ((ulong) 1 << 21) + ((ulong) 1 << 42);
   exponents3[3] = (ulong) 2 + ((ulong) 2 << 21) + ((ulong) 1 << 42);
   exponents3[4] = (ulong) 3 + ((ulong) 2 << 21) + ((ulong) 1 << 42);
   exponents3[5] = (ulong) 0 + ((ulong) 0 << 21) + ((ulong) 2 << 42);
   exponents3[6] = (ulong) 1 + ((ulong) 0 << 21) + ((ulong) 2 << 42);
   exponents3[7] = (ulong) 2 + ((ulong) 0 << 21) + ((ulong) 2 << 42);
   exponents3[8] = (ulong) 1 + ((ulong) 1 << 21) + ((ulong) 2 << 42);
   exponents3[9] = (ulong) 3 + ((ulong) 1 << 21) + ((ulong) 2 << 42);
   exponents3[10] = (ulong) 2 + ((ulong) 2 << 21) + ((ulong) 2 << 42);
   exponents3[11] = (ulong) 3 + ((ulong) 2 << 21) + ((ulong) 2 << 42);
   exponents3[12] = (ulong) 4 + ((ulong) 2 << 21) + ((ulong) 2 << 42);
   exponents3[13] = (ulong) 2 + ((ulong) 1 << 21) + ((ulong) 3 << 42);

   ulong monomial3 = (ulong) 2 + ((ulong) 1 << 21) + ((ulong) 2 << 42);
   
   get_period_sequence(zeroCoefficients3, coefficients3, exponents3, 14, monomial3, 10, primes3, 10);

   for(i=0; i<100; i++){
      fmpz_clear(zeroCoefficients3[i]);
   }
*/
 
   printf("PASS\n");
   return 0;
}
