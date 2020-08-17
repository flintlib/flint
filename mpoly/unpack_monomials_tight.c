/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"

void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len,
                                          slong * mults, slong num, slong bits)
{
   slong i, j;
   ulong exp;
   slong * prods;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 1; i <= num; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   for (i = 0; i < len; i++)
   {
      exp = 0;
         
      for (j = 0; j < num; j++)
         exp += (e2[i] % prods[j + 1])/prods[j] << bits*j;

      e1[i] = exp;

   }

   TMP_END;
}

