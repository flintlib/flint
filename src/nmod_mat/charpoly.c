/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

void
_nmod_mat_charpoly_berkowitz(nn_ptr cp, const nmod_mat_t mat, nmod_t mod)
{
    const slong n = mat->r;

    if (mod.n == 1)
    {
        _nmod_vec_zero(cp, n + 1);
    }
    else if (n == 0)
    {
        cp[0] = 1;
    }
    else if (n == 1)
    {
        cp[0] = nmod_neg(nmod_mat_entry(mat, 0, 0), mod);
        cp[1] = 1;
    }
    else if (n == 2)
    {
        cp[0] = nmod_sub(nmod_mul(nmod_mat_entry(mat, 0, 0), nmod_mat_entry(mat, 1, 1), mod),
                         nmod_mul(nmod_mat_entry(mat, 0, 1), nmod_mat_entry(mat, 1, 0), mod), mod);
        cp[1] = nmod_add(nmod_mat_entry(mat, 0, 0), nmod_mat_entry(mat, 1, 1), mod);
        cp[1] = nmod_neg(cp[1], mod);
        cp[2] = 1;
    }
    else
    {
        slong i, k, t;
        nn_ptr a, A, s;
        TMP_INIT;

        TMP_START;
        a = TMP_ALLOC(sizeof(ulong) * (n * n));
        A = a + (n - 1) * n;

        const dot_params_t params = _nmod_vec_dot_params(n, mod);

        _nmod_vec_zero(cp, n + 1);
        cp[0] = nmod_neg(nmod_mat_entry(mat, 0, 0), mod);

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                a[0 * n + i] = nmod_mat_entry(mat, i, t);
            }

            A[0] = nmod_mat_entry(mat, t, t);

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    s[0] = _nmod_vec_dot(mat->rows[i], a + (k - 1) * n, t + 1, mod, params);
                }

                A[k] = a[k * n + t];
            }

            A[t] = _nmod_vec_dot(mat->rows[t], a + (t - 1) * n, t + 1, mod, params);

            for (k = 0; k <= t; k++)
            {
                cp[k] = nmod_sub(cp[k], _nmod_vec_dot_rev(A, cp, k, mod, params), mod);
                cp[k] = nmod_sub(cp[k], A[k], mod);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            cp[i] = cp[i - 1];

        cp[0] = 1;
        _nmod_poly_reverse(cp, cp, n + 1, n + 1);

        TMP_END;
    }
}

void nmod_mat_charpoly_berkowitz(nmod_poly_t cp, const nmod_mat_t mat)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_mat_charpoly_berkowitz).  Non-square matrix.\n");
    }

    nmod_poly_fit_length(cp, mat->r + 1);
    _nmod_poly_set_length(cp, mat->r + 1);
    _nmod_mat_charpoly_berkowitz(cp->coeffs, mat, mat->mod);
}

void nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M)
{
   slong n = M->r, i, j, k;
   ulong ** A;
   ulong * V, * W, * T;
   ulong h;
   nmod_poly_t b;
   nmod_mat_t M2;
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
   const dot_params_t params = _nmod_vec_dot_params(n, p->mod);
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

         A[n - i - 1][j - 1] = _nmod_vec_dot(T, W, n - i, p->mod, params);
      }

      for (j = n - i; j <= n - 1; j++)
      {
         for (k = 1; k <= n - i; k++)
            T[k - 1] = A[k - 1][j - 1];

         A[n - i - 1][j - 1] = n_addmod(_nmod_vec_dot(T, W, n - i, p->mod, params), W[j], p->mod.n);
      }

      for (k = 1; k <= n - i; k++)
         T[k - 1] = A[k - 1][j - 1];

      A[n - i - 1][n - 1] = _nmod_vec_dot(T, W, n - i, p->mod, params);

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

void
nmod_mat_charpoly(nmod_poly_t cp, const nmod_mat_t mat)
{
    if (mat->r <= 8 || !n_is_prime(mat->mod.n))
        nmod_mat_charpoly_berkowitz(cp, mat);
    else
        nmod_mat_charpoly_danilevsky(cp, mat);
}
