/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"

int main(void)
{
   int i;
   FLINT_TEST_INIT(state);
   
   flint_printf("ll_knuth_schroeppel....");
   fflush(stdout);
 
   

   for (i = 0; i < 10000; i++) /* Test random n */
   {
      mp_limb_t hi = 0, lo;
      qs_t qs_inf;
      mp_bitcnt_t bits;
      
      bits = n_randint(state, 2*FLINT_BITS) + 1;
      if (bits > FLINT_BITS)
      {
          lo = n_randlimb(state);
          hi = n_randbits(state, bits - FLINT_BITS);
      } else
          lo = n_randbits(state, bits);
      
      qsieve_ll_init(qs_inf, hi, lo);
      qsieve_ll_knuth_schroeppel(qs_inf);
      qsieve_ll_clear(qs_inf);
   }
   
   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
