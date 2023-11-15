/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "qfb.h"

slong qfb_hash_find(qfb_hash_t * qhash, qfb_t q, slong depth)
{
   slong size = (1L<<depth), i;
   fmpz_t r;

   fmpz_init(r);

   fmpz_fdiv_r_2exp(r, q->a, depth);
   i = fmpz_get_ui(r);

   while (!fmpz_is_zero(qhash[i].q->a))
   {
      if (fmpz_cmp(qhash[i].q->a, q->a) == 0)
      {
         if (fmpz_cmpabs(qhash[i].q->b, q->b) == 0)
         {
            fmpz_clear(r);
            return i;
         }
      }

      i++;
      if (i == size)
         i = 0;
   }

   fmpz_clear(r);
   return -1;
}
