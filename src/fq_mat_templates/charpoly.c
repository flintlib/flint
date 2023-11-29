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
TEMPLATE(T, mat_charpoly_danilevsky) (TEMPLATE(T, poly_t) p,
                      const TEMPLATE(T, mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
   slong n = M->r, i, j, k;
   TEMPLATE(T, struct) * V, * W, * T;
   TEMPLATE(T, poly_t) b;
   TEMPLATE(T, t) c, h;
   TEMPLATE(T, mat_t) A;

   if (n == 0)
   {
      TEMPLATE(T, poly_one) (p, ctx);
      return;
   }

   TEMPLATE(T, init) (c, ctx);

   if (n == 1)
   {
      TEMPLATE(T, set_ui) (c, 1, ctx);
      TEMPLATE(T, poly_set_coeff) (p, 1, c, ctx);
      TEMPLATE(T, neg) (c, TEMPLATE(T, mat_entry) (M, 0, 0), ctx);
      TEMPLATE(T, poly_set_coeff) (p, 0, c, ctx);
      _TEMPLATE(T, poly_set_length) (p, 2, ctx);
      TEMPLATE(T, clear) (c, ctx);
      return;
   }

   i = 1;
   TEMPLATE(T, init) (h, ctx);
   TEMPLATE(T, poly_one) (p, ctx);
   TEMPLATE(T, poly_init) (b, ctx);
   TEMPLATE(T, mat_init_set) (A, M, ctx);
   V = (TEMPLATE(T, struct) *) _TEMPLATE(T, vec_init) (n, ctx);
   W = (TEMPLATE(T, struct) *) _TEMPLATE(T, vec_init) (n, ctx);
   T = (TEMPLATE(T, struct) *) _TEMPLATE(T, vec_init) (n, ctx);

   while (i < n)
   {
      TEMPLATE(T, set) (h, TEMPLATE(T, mat_entry) (A, n - i, n - i - 1), ctx);

      while (TEMPLATE(T, is_zero) (h, ctx))
      {
         k = 1;
         while (k < n - i && TEMPLATE(T, is_zero) (TEMPLATE(T, mat_entry) (A, n - i, n - i - k - 1), ctx))
            k++;

         if (k == n - i)
         {
            TEMPLATE(T, poly_fit_length) (b, i + 1, ctx);
            TEMPLATE(T, set_ui) (c, 1, ctx);
            TEMPLATE(T, poly_set_coeff) (b, i, c, ctx);
            for (k = 1; k <= i; k++)
            {
               TEMPLATE(T, neg) (c, TEMPLATE(T, mat_entry) (A, n - i, n - k), ctx);
               TEMPLATE(T, poly_set_coeff) (b, k - 1, c, ctx);
            }
            _TEMPLATE(T, poly_set_length) (b, i + 1, ctx);
            TEMPLATE(T, poly_mul) (p, p, b, ctx);

            n -= i;
            i = 1;

            if (n == 1)
            {
               TEMPLATE(T, set_ui) (c, 1, ctx);
               TEMPLATE(T, poly_set_coeff) (b, 1, c, ctx);
               TEMPLATE(T, neg) (c, TEMPLATE(T, mat_entry) (A, 0, 0), ctx);
               TEMPLATE(T, poly_set_coeff) (b, 0, c, ctx);
               _TEMPLATE(T, poly_set_length) (b, 2, ctx);
               TEMPLATE(T, poly_mul) (p, p, b, ctx);

               goto cleanup;
            }
         } else
         {
            TEMPLATE(T, struct) * ptr;

            ptr = A->rows[n - i - k - 1];
            A->rows[n - i - k - 1] = A->rows[n - i - 1];
            A->rows[n - i - 1] = ptr;

            for (j = 1; j <= n - i + 1; j++)
            {
               TEMPLATE(T, swap) (TEMPLATE(T, mat_entry) (A, j - 1, n - i - k - 1),
                                  TEMPLATE(T, mat_entry) (A, j - 1, n - i - 1), ctx);
            }
         }

         TEMPLATE(T, set) (h, TEMPLATE(T, mat_entry) (A, n - i, n - i - 1), ctx);
      }

      TEMPLATE(T, neg) (h, h, ctx);
      TEMPLATE(T, inv) (h, h, ctx);

      for (j = 1; j <= n; j++)
      {
         TEMPLATE(T, mul) (V + j - 1, TEMPLATE(T, mat_entry) (A, n - i, j - 1), h, ctx);
         TEMPLATE(T, set) (W + j - 1, TEMPLATE(T, mat_entry) (A, n - i, j - 1), ctx);
      }

      TEMPLATE(T, neg) (h, h, ctx);

      for (j = 1; j <= n - i; j++)
      {
         for (k = 1; k <= n - i - 1; k++)
         {
            TEMPLATE(T, mul) (c, TEMPLATE(T, mat_entry) (A, j - 1, n - i - 1), V + k - 1, ctx);
            TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, j - 1, k - 1),
                              TEMPLATE(T, mat_entry) (A, j - 1, k - 1), c, ctx);
         }

         for (k = n - i + 1; k <= n; k++)
         {
            TEMPLATE(T, mul) (c, TEMPLATE(T, mat_entry) (A, j - 1, n - i - 1), V + k - 1, ctx);
            TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, j - 1, k - 1),
                              TEMPLATE(T, mat_entry) (A, j - 1, k - 1), c, ctx);
         }

         TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (A, j - 1, n - i - 1),
                           TEMPLATE(T, mat_entry) (A, j - 1, n - i - 1), h, ctx);
      }

      for (j = 1; j <= n - i - 1; j++)
      {
         TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1),
                           TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1), W + n - i - 1, ctx);

         for (k = 1; k < n - i; k++)
         {
            TEMPLATE(T, mul) (c, TEMPLATE(T, mat_entry) (A, k - 1, j - 1), W + k - 1, ctx);
            TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1),
                              TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1), c, ctx);
         }
      }

      for (j = n - i; j <= n - 1; j++)
      {
         TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1),
                           TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1), W + n - i - 1, ctx);

         for (k = 1; k < n - i; k++)
         {
            TEMPLATE(T, mul) (c, TEMPLATE(T, mat_entry) (A, k - 1, j - 1), W + k - 1, ctx);
            TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1),
                              TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1), c, ctx);
         }

         TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1),
                              TEMPLATE(T, mat_entry) (A, n - i - 1, j - 1), W + j, ctx);
      }

      TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (A, n - i - 1, n - 1),
                        TEMPLATE(T, mat_entry) (A, n - i - 1, n - 1), W + n - i - 1, ctx);

      for (k = 1; k < n - i; k++)
      {
         TEMPLATE(T, mul) (c, TEMPLATE(T, mat_entry) (A, k - 1, n - 1), W + k - 1, ctx);
         TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (A, n - i - 1, n - 1),
                           TEMPLATE(T, mat_entry) (A, n - i - 1, n - 1), c, ctx);
      }

      i++;
   }

   TEMPLATE(T, poly_fit_length) (b, n + 1, ctx);
   TEMPLATE(T, set_ui) (c, 1, ctx);
   TEMPLATE(T, poly_set_coeff) (b, n, c, ctx);
   for (i = 1; i <= n; i++)
   {
      TEMPLATE(T, neg) (c, TEMPLATE(T, mat_entry) (A, 0, n - i), ctx);
      TEMPLATE(T, poly_set_coeff) (b, i - 1, c, ctx);
   }
   _TEMPLATE(T, poly_set_length) (b, n + 1, ctx);
   TEMPLATE(T, poly_mul) (p, p, b, ctx);

cleanup:

   TEMPLATE(T, mat_clear) (A, ctx);
   TEMPLATE(T, clear) (c, ctx);
   TEMPLATE(T, clear) (h, ctx);
   TEMPLATE(T, poly_clear) (b, ctx);
   _TEMPLATE(T, vec_clear) (T, A->r, ctx);
   _TEMPLATE(T, vec_clear) (V, A->r, ctx);
   _TEMPLATE(T, vec_clear) (W, A->r, ctx);
}

void
TEMPLATE(T, mat_charpoly)(TEMPLATE(T, poly_t) p,
                       const TEMPLATE(T, mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
   TEMPLATE(T, mat_t) A;

   TEMPLATE(T, mat_init) (A, M->r, M->c, ctx);
   TEMPLATE(T, mat_set) (A, M, ctx);

   if (A->r != A->c)
   {
       flint_throw(FLINT_ERROR, "Exception (fq_mat_charpoly).  Non-square matrix.\n");
   }

   TEMPLATE(T, mat_charpoly_danilevsky) (p, A, ctx);

   TEMPLATE(T, mat_clear) (A, ctx);
}

#endif
