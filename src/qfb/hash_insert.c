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

void qfb_hash_insert(qfb_hash_t * qhash, qfb_t q, qfb_t q2, slong iter, slong depth)
{
   slong size = (1L<<depth), i;
   fmpz_t r;

   fmpz_init(r);

   fmpz_fdiv_r_2exp(r, q->a, depth);
   i = fmpz_get_ui(r);

   while (!fmpz_is_zero(qhash[i].q->a))
   {
      i++;
      if (i == size)
         i = 0;
   }

   qfb_set(qhash[i].q, q);
   qhash[i].iter = iter;
   if (q2 != NULL)
      qfb_set(qhash[i].q2, q2);

   fmpz_clear(r);
}
