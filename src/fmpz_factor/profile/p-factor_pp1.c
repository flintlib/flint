/*
    Copyright 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"

int main(void)
{
   fmpz_t n, p;
   ulong c;
   ulong B1;

   FLINT_TEST_INIT(state);

   fmpz_init(n);
   fmpz_init(p);

   while(1)
   {
      flint_printf("Enter number to be factored: "); fflush(stdout);
      if (!fmpz_read(n))
      {
         flint_printf("Read failed\n");
         flint_abort();
      }

      flint_printf("Enter B1: "); fflush(stdout);
      if (!flint_scanf("%wu", &B1))
      {
         flint_printf("Read failed\n");
         flint_abort();
      }

      do
      {
         c = n_randlimb(state);
      } while (c <= UWORD(2));

      if (fmpz_factor_pp1(p, n, B1, B1/100, c))
      {
         flint_printf("Factor: ");
         fmpz_print(p);
         flint_printf("\n");
      } else
         flint_printf("Factor not found!\n");
   } while(1);

   flint_rand_clear(state);

   fmpz_clear(n);
   fmpz_clear(p);

   return 0;
}
