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

void TEMPLATE(T, mat_similarity) (TEMPLATE(T, mat_t) A, slong r,
                                TEMPLATE(T, t) d, const TEMPLATE(T, ctx_t) ctx)
{
   slong n = A->r, i, j;
   TEMPLATE(T, t) t;

   TEMPLATE(T, init) (t, ctx);

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
      {
         TEMPLATE(T, mul) (t, TEMPLATE(T, mat_entry) (A, i, r), d, ctx);
         TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, i, j), TEMPLATE(T, mat_entry) (A, i, j), t, ctx);
      }

      for (j = r + 1; j < n; j++)
      {
         TEMPLATE(T, mul) (t, TEMPLATE(T, mat_entry) (A, i, r), d, ctx);
         TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, i, j), TEMPLATE(T, mat_entry) (A, i, j), t, ctx);
      }
   }

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
      {
         TEMPLATE(T, mul) (t, TEMPLATE(T, mat_entry) (A, j, i), d, ctx);
         TEMPLATE(T, sub) (TEMPLATE(T, mat_entry) (A, r, i), TEMPLATE(T, mat_entry) (A, r, i), t, ctx);
      }

      for (j = r + 1; j < n; j++)
      {
         TEMPLATE(T, mul) (t, TEMPLATE(T, mat_entry) (A, j, i), d, ctx);
         TEMPLATE(T, sub) (TEMPLATE(T, mat_entry) (A, r, i), TEMPLATE(T, mat_entry) (A, r, i), t, ctx);
      }
   }

   TEMPLATE(T, clear) (t, ctx);
}

#endif
