/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

void nmod_mat_minpoly_with_gens(nmod_poly_t p, const nmod_mat_t X, ulong * P)
{
   slong n = X->r, i, j, c, c1, c2, r1, r2;
   ulong ** A, ** B, ** v, t, h;
   slong  * P1, * P2, * L1, * L2;
   nmod_mat_t matA, matB, matv;
   int first_poly = 1, indep = 1;
   nmod_poly_t b, g;
   TMP_INIT;

   if (X->r != X->c)
   {
       flint_throw(FLINT_ERROR, "Exception (nmod_mat_charpoly).  Non-square matrix.\n");
   }

   if (n == 0)
   {
      nmod_poly_one(p);
      return;
   }

   if (n == 1)
   {
      nmod_poly_set_coeff_ui(p, 1, 1);
      nmod_poly_set_coeff_ui(p, 0, n_negmod(X->rows[0][0], p->mod.n));
      _nmod_poly_set_length(p, 2);
      if (P != NULL)
         P[0] = 1;
      return;
   }

   TMP_START;

   nmod_poly_init(b, p->mod.n);
   nmod_poly_init(g, p->mod.n);
   nmod_poly_one(p);
   nmod_mat_init(matA, n + 1, 2*n + 1, p->mod.n);
   nmod_mat_init(matB, n, n, p->mod.n);
   nmod_mat_init(matv, n, 1, p->mod.n);

   A = matA->rows;
   B = matB->rows;
   v = matv->rows;

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
         v[i][0] = 0;
         B[r2][i] = 0;
         A[0][i] = 0;
      }

      P1[c2] = 0;
      P2[c2] = r2;

      v[c2][0] = 1;
      B[r2][c2] = 1;
      A[0][c2] = 1;
      A[0][n] = 1;
      if (P != NULL)
         P[c2] = 1;

      indep = 1;

      r1 = 0;
      c1 = -WORD(1);

      while (c1 < n && r1 < n)
      {
         r1++;
         r2 = indep ? r2 + 1 : r2;

         nmod_mat_mul(matv, X, matv);
         v = matv->rows;

         for (i = 0; i < n; i++)
            A[r1][i] = v[i][0];

         for (i = n; i < n + r1; i++)
            A[r1][i] = 0;

         A[r1][n + r1] = 1;

         c1 = nmod_mat_reduce_row(matA, P1, L1, r1);

         if (indep && r2 < n && !first_poly)
         {
            for (i = 0; i < n; i++)
               B[r2][i] = v[i][0];

            c = nmod_mat_reduce_row(matB, P2, L2, r2);

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

      nmod_poly_fit_length(b, r1 + 1);

      h = n_invmod(A[r1][n + r1], p->mod.n);

      for (i = 0; i < r1 + 1; i++)
      {
         t = n_mulmod2_preinv(A[r1][n + i], h, p->mod.n, p->mod.ninv);
         nmod_poly_set_coeff_ui(b, i, t);
      }
      _nmod_poly_set_length(b, r1 + 1);

      nmod_poly_gcd(g, p, b);
      nmod_poly_mul(p, p, b);
      nmod_poly_div(p, p, g);

      if (first_poly && r2 < n)
      {
         for (i = 0; i < r1; i++)
         {
            for (j = 0; j < n; j++)
               B[i][j] = A[i][j];
         }
      }

      first_poly = 0;
   }

   nmod_mat_clear(matA);
   nmod_mat_clear(matB);
   nmod_mat_clear(matv);

   nmod_poly_clear(b);
   nmod_poly_clear(g);

   TMP_END;
}

void nmod_mat_minpoly(nmod_poly_t p, const nmod_mat_t X)
{
   nmod_mat_minpoly_with_gens(p, X, NULL);
}
