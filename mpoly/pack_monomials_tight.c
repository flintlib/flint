/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"

void mpoly_pack_monomials_tight(ulong * exp1, const ulong * exp2,
            slong len, const slong * mults, slong num, slong extra, slong bits)
{
   slong i, j;
   ulong e1, e2;
   ulong mask = (UWORD(1) << bits) - 1;
   slong shift = FLINT_BITS - (num + extra)*bits;

   for (i = 0; i < len; i++)
   {
      e2 = exp2[i] >> shift;
      e1 = ((e2 >> (num - 1)*bits) & mask);
      
      for (j = num - 2; j >= 0; j--)
      {
         e1 *= mults[j];
         e1 += ((e2 >> j*bits) & mask);
      }

      exp1[i] = e1;
   }
}


