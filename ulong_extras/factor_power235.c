/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2009 Thomas Boothby

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with standard libraries */
#include <math.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_factor_power235(ulong * exp, mp_limb_t n)
{
   static char mod63[63] = {7,7,4,0,5,4,0,5,6,5,4,4,0,4,4,0,5,4,5,4,4,0,
                            5,4,0,5,4,6,7,4,0,4,4,0,4,6,7,5,4,0,4,4,0,5,
                            4,4,5,4,0,5,4,0,4,4,4,6,4,0,5,4,0,4,6};
   static char mod61[61] = {7,7,0,3,1,1,0,0,2,3,0,6,1,5,5,1,1,0,0,1,3,4,
                            1,2,2,1,0,3,2,4,0,0,4,2,3,0,1,2,2,1,4,3,1,0,
                            0,1,1,5,5,1,6,0,3,2,0,0,1,1,3,0,7};
   static char mod44[44] = {7,7,0,2,3,3,0,2,2,3,0,6,7,2,0,2,3,2,0,2,3,6,
                            0,6,2,3,0,2,2,2,0,2,6,7,0,2,3,3,0,2,2,2,0,6};
   static char mod31[31] = {7,7,3,0,3,5,4,1,3,1,1,0,0,0,1,2,3,0,1,1,1,0,
                            0,2,0,5,4,2,1,2,6};
   char t;
   
   t = mod31[n%31];
   if (!t) return UWORD(0);

   t &= mod44[n%44];
   if (!t) return UWORD(0);

   t &= mod61[n%61];
   if (!t) return UWORD(0);

   t&= mod63[n%63];

   if (t & 1) 
   {
      double x = sqrt((double) n);
      mp_limb_t y = (mp_limb_t) (x + 0.5);
      if (n == n_pow(y, 2))
      {
         *exp = 2;
         return y;
      }
   }
    
   if (t & 2) 
   {
      double x = pow((double) n, 1.0 / 3.0);
      mp_limb_t y = (mp_limb_t) (x + 0.5);
      if (n == n_pow(y, 3))
      {
         *exp = 3;
         return y;
      }
   }
    
   if (t & 4) 
   {
      double x = pow((double) n, 1.0 / 5.0);
      mp_limb_t y = (mp_limb_t) (x + 0.5);
      if (n == n_pow(y, 5))
      {
         *exp = 5;
         return y;
      }
   }

   return UWORD(0);
}
