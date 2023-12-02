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

void
TEMPLATE(T, mat_minpoly) (TEMPLATE(T, poly_t) p,
                      const TEMPLATE(T, mat_t) X, const TEMPLATE(T, ctx_t) ctx)
{
   slong n = X->r, i, j, c, c1, c2, r1, r2;
   slong  * P1, * P2, * L1, * L2;
   TEMPLATE(T, mat_t) A, B, v;
   int first_poly = 1, indep = 1;
   TEMPLATE(T, poly_t) b, g, r;
   TEMPLATE(T, t) t, h;
   TMP_INIT;

   if (X->r != X->c)
   {
       flint_throw(FLINT_ERROR, "Exception (fq_mat_charpoly).  Non-square matrix.\n");
   }

   if (n == 0)
   {
      TEMPLATE(T, poly_one) (p, ctx);
      return;
   }

   TEMPLATE(T, init) (t, ctx);

   if (n == 1)
   {
      TEMPLATE(T, set_ui) (t, 1, ctx);
      TEMPLATE(T, poly_set_coeff) (p, 1, t, ctx);
      TEMPLATE(T, neg) (t, TEMPLATE(T, mat_entry) (X, 0, 0), ctx);
      TEMPLATE(T, poly_set_coeff) (p, 0, t, ctx);
      _TEMPLATE(T, poly_set_length) (p, 2, ctx);
      TEMPLATE(T, clear) (t, ctx);
      return;
   }

   TMP_START;

   TEMPLATE(T, init) (h, ctx);
   TEMPLATE(T, poly_init) (b, ctx);
   TEMPLATE(T, poly_init) (g, ctx);
   TEMPLATE(T, poly_init) (r, ctx);
   TEMPLATE(T, poly_one) (p, ctx);
   TEMPLATE(T, mat_init) (A, n + 1, 2*n + 1, ctx);
   TEMPLATE(T, mat_init) (B, n, n, ctx);
   TEMPLATE(T, mat_init) (v, n, 1, ctx);

   L1 = (slong *) TMP_ALLOC((n + 1)*sizeof(slong));
   L2 = (slong *) TMP_ALLOC(n*sizeof(slong));
   P1 = (slong *) TMP_ALLOC((2*n + 1)*sizeof(slong));
   P2 = (slong *) TMP_ALLOC(n*sizeof(slong));

   for (i = 1; i <= n + 1; i++)
      L1[i - 1] = n + i;

   for (i = 1; i <= n; i++)
      L2[i - 1] = n;

   for (i = 1; i < n; i++)
      P2[i] = -WORD(1);
   P2[0] = 0;

   r2 = c2 = 0;
   first_poly = 1;

   while (r2 < n)
   {
      for (i = 0; i < 2*n + 1; i++)
         P1[i] = -WORD(1);

      for (i = 0; i < n; i++)
      {
         TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (v, i, 0), ctx);
         TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (B, r2, i), ctx);
         TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (A, 0, i), ctx);
      }

      P1[c2] = 0;
      P2[c2] = r2;

      TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (v, c2, 0), ctx);
      TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (B, r2, c2), ctx);
      TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (A, 0, c2), ctx);
      TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (A, 0, n), ctx);

      indep = 1;

      r1 = 0;
      c1 = -WORD(1);

      while (c1 < n && r1 < n)
      {
         r1++;
         r2 = indep ? r2 + 1 : r2;

         TEMPLATE(T, mat_mul) (v, X, v, ctx);

         for (i = 0; i < n; i++)
            TEMPLATE(T, set) (TEMPLATE(T, mat_entry) (A, r1, i), TEMPLATE(T, mat_entry) (v, i, 0), ctx);

         for (i = n; i < n + r1; i++)
            TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (A, r1, i), ctx);

         TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (A, r1, n + r1), ctx);

         c1 = TEMPLATE(T, mat_reduce_row) (A, P1, L1, r1, ctx);

         if (indep && r2 < n && !first_poly)
         {
            for (i = 0; i < n; i++)
               TEMPLATE(T, set) (TEMPLATE(T, mat_entry) (B, r2, i), TEMPLATE(T, mat_entry) (v, i, 0), ctx);

            c = TEMPLATE(T, mat_reduce_row) (B, P2, L2, r2, ctx);

            indep = c != -WORD(1);
         }
      }

      if (first_poly)
      {
         for (i = 0; i < n; i++)
            P2[i] = P1[i];

         r2 = r1;
      }

      c = -WORD(1);

      for (i = c2 + 1; i < n; i++)
      {
         if (P2[i] == -WORD(1))
         {
            c = i;
            break;
         }
      }

      c2 = c;

      TEMPLATE(T, poly_fit_length) (b, r1 + 1, ctx);

      TEMPLATE(T, inv) (h, TEMPLATE(T, mat_entry) (A, r1, n + r1), ctx);

      for (i = 0; i < r1 + 1; i++)
      {
         TEMPLATE(T, mul) (t, TEMPLATE(T, mat_entry) (A, r1, n + i), h, ctx);
         TEMPLATE(T, poly_set_coeff) (b, i, t, ctx);
      }
      _TEMPLATE(T, poly_set_length) (b, r1 + 1, ctx);

      TEMPLATE(T, poly_gcd) (g, p, b, ctx);
      TEMPLATE(T, poly_mul) (p, p, b, ctx);
      TEMPLATE(T, poly_divrem) (p, r, p, g, ctx);

      if (first_poly && r2 < n)
      {
         for (i = 0; i < r1; i++)
         {
            for (j = 0; j < n; j++)
               TEMPLATE(T, set) (TEMPLATE(T, mat_entry) (B, i, j),  TEMPLATE(T, mat_entry) (A, i, j), ctx);
         }
      }

      first_poly = 0;
   }

   TEMPLATE(T, mat_clear) (A, ctx);
   TEMPLATE(T, mat_clear) (B, ctx);
   TEMPLATE(T, mat_clear) (v, ctx);

   TEMPLATE(T, poly_clear) (b, ctx);
   TEMPLATE(T, poly_clear) (g, ctx);
   TEMPLATE(T, poly_clear) (r, ctx);

   TEMPLATE(T, clear) (t, ctx);
   TEMPLATE(T, clear) (h, ctx);

   TMP_END;
}

#endif
