/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#include "templates.h"

slong TEMPLATE(T, mat_reduce_row)(TEMPLATE(T, mat_t) A, slong * P, slong * L, 
                                         slong m, const TEMPLATE(T, ctx_t) ctx)
{
   slong n = A->c, i, j, r;
   TEMPLATE(T, t) h;

   TEMPLATE(T, init) (h, ctx);

   for (i = 0; i < n; i++)
   {
      if (!TEMPLATE(T, is_zero) (TEMPLATE(T, mat_entry) (A, m, i), ctx))
      {
         r = P[i];
         if (r != -WORD(1))
         {
            for (j = i + 1; j < L[r]; j++)
            {
               TEMPLATE(T, mul) (h, TEMPLATE(T, mat_entry) (A, r, j), TEMPLATE(T, mat_entry) (A, m, i), ctx);
               TEMPLATE(T, sub) (TEMPLATE(T, mat_entry) (A, m, j), TEMPLATE(T, mat_entry) (A, m, j), h, ctx);
            }
 
            TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (A, m, i), ctx);
         } else
         {
            TEMPLATE(T, inv) (h, TEMPLATE(T, mat_entry) (A, m, i), ctx);
            TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (A, m, i), ctx);
           
            for (j = i + 1; j < L[m]; j++)
               TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (A, m, j), TEMPLATE(T, mat_entry) (A, m, j), h, ctx);
               
            P[i] = m;

            TEMPLATE(T, clear) (h, ctx);
       
            return i;
         }
      }
   }

   TEMPLATE(T, clear) (h, ctx);
   
   return -WORD(1);
}

#endif
