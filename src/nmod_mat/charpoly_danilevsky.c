/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

void nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M)
{
   slong n = M->r, i, j, k;
   ulong ** A;
   ulong * V, * W, * T;
   ulong h;
   nmod_poly_t b;
   nmod_mat_t M2;
   int num_limbs;
   TMP_INIT;

   if (M->r != M->c)
   {
      flint_throw(FLINT_ERROR, "Exception (nmod_mat_charpoly_danilevsky).  Non-square matrix.\n");
   }

   if (n == 0)
   {
      nmod_poly_one(p);
      return;
   }

   if (n == 1)
   {
      nmod_poly_set_coeff_ui(p, 1, 1);
      nmod_poly_set_coeff_ui(p, 0, n_negmod(M->rows[0][0], p->mod.n));
      _nmod_poly_set_length(p, 2);
      return;
   }

   TMP_START;

   i = 1;
   num_limbs = _nmod_vec_dot_bound_limbs(n, p->mod);
   nmod_poly_one(p);
   nmod_poly_init(b, p->mod.n);
   nmod_mat_init_set(M2, M);
   V = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   W = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   T = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   A = M2->rows;

   while (i < n)
   {
      h = A[n - i][n - i - 1];

      while (h == 0)
      {
         k = 1;
         while (k < n - i && A[n - i][n - i - k - 1] == 0)
            k++;

         if (k == n - i)
         {
            nmod_poly_fit_length(b, i + 1);
            nmod_poly_set_coeff_ui(b, i, 1);
            for (k = 1; k <= i; k++)
               nmod_poly_set_coeff_ui(b, k - 1, n_negmod(A[n - i][n - k], p->mod.n));
            _nmod_poly_set_length(b, i + 1);
            nmod_poly_mul(p, p, b);
            n -= i;
            i = 1;

            if (n == 1)
            {
               nmod_poly_set_coeff_ui(b, 1, 1);
               nmod_poly_set_coeff_ui(b, 0, n_negmod(A[0][0], p->mod.n));
               _nmod_poly_set_length(b, 2);
               nmod_poly_mul(p, p, b);
               goto cleanup;
            }
         } else
         {
            ulong * ptr;
            ulong t;

            ptr = A[n - i - k - 1];
            A[n - i - k - 1] = A[n - i - 1];
            A[n - i - 1] = ptr;

            for (j = 1; j <= n - i + 1; j++)
            {
               t = A[j - 1][n - i - k - 1];
               A[j - 1][n - i - k - 1] = A[j - 1][n - i - 1];
               A[j - 1][n - i - 1] = t;
            }
         }

         h = A[n - i][n - i - 1];
      }

      h = n_invmod(n_negmod(h, p->mod.n), p->mod.n);

      for (j = 1; j <= n; j++)
      {
         V[j - 1] = n_mulmod2_preinv(A[n - i][j - 1], h, p->mod.n, p->mod.ninv);
         W[j - 1] = A[n - i][j - 1];
      }

      h = n_negmod(h, p->mod.n);

      for (j = 1; j <= n - i; j++)
      {
         for (k = 1; k <= n - i - 1; k++)
            NMOD_ADDMUL(A[j - 1][k - 1], A[j - 1][n - i - 1], V[k - 1], p->mod);

         for (k = n - i + 1; k <= n; k++)
            NMOD_ADDMUL(A[j - 1][k - 1], A[j - 1][n - i - 1], V[k - 1], p->mod);

         A[j - 1][n - i - 1] = n_mulmod2_preinv(A[j - 1][n - i - 1], h, p->mod.n, p->mod.ninv);
      }

      for (j = 1; j <= n - i - 1; j++)
      {
         for (k = 1; k <= n - i; k++)
            T[k - 1] = A[k - 1][j - 1];

         A[n - i - 1][j - 1] = _nmod_vec_dot(T, W, n - i, p->mod, num_limbs);
      }

      for (j = n - i; j <= n - 1; j++)
      {
         for (k = 1; k <= n - i; k++)
            T[k - 1] = A[k - 1][j - 1];

         A[n - i - 1][j - 1] = n_addmod(_nmod_vec_dot(T, W, n - i, p->mod, num_limbs), W[j], p->mod.n);
      }

      for (k = 1; k <= n - i; k++)
         T[k - 1] = A[k - 1][j - 1];

      A[n - i - 1][n - 1] = _nmod_vec_dot(T, W, n - i, p->mod, num_limbs);

      i++;
   }

   nmod_poly_fit_length(b, n + 1);
   nmod_poly_set_coeff_ui(b, n, 1);
   for (i = 1; i <= n; i++)
      nmod_poly_set_coeff_ui(b, i - 1, n_negmod(A[0][n - i], p->mod.n));
   _nmod_poly_set_length(b, n + 1);
   nmod_poly_mul(p, p, b);

cleanup:

   nmod_mat_clear(M2);
   nmod_poly_clear(b);
   TMP_END;
}

