/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "qfb.h"

void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D)
{
   int done = 0;
   fmpz_t t;

   qfb_set(r, f);
   
   fmpz_init(t);

   while(!done)
   {
      done = 1;

      if (fmpz_cmp(r->c, r->a) < 0)
      {
         fmpz_swap(r->a, r->c);
         fmpz_neg(r->b, r->b);

         done = 0;
      }

      if (fmpz_cmpabs(r->b, r->a) > 0)
      {
         fmpz_add(t, r->a, r->a);
         fmpz_fdiv_r(r->b, r->b, t);
         if (fmpz_cmp(r->b, r->a) > 0)
            fmpz_sub(r->b, r->b, t);

         fmpz_add(t, t, t);
         fmpz_mul(r->c, r->b, r->b);
         fmpz_sub(r->c, r->c, D);
         fmpz_divexact(r->c, r->c, t);

         done = 0;
      }
   }

   if (fmpz_cmpabs(r->a, r->b) == 0 || fmpz_cmp(r->a, r->c) == 0)
      if (fmpz_sgn(r->b) < 0)
         fmpz_neg(r->b, r->b);

   fmpz_clear(t);
}
