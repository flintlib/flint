/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

slong nmod_mat_reduce_row(nmod_mat_t M, slong * P, slong * L, slong m)
{
   slong n = M->c, i, j, r, bits, limbs;
   ulong ** A = M->rows;
   ulong h, hi, lo;
   ulong * rowm;
   slong res = -WORD(1);

   TMP_INIT;

   bits = FLINT_BIT_COUNT(M->mod.n)*2 + FLINT_BIT_COUNT(m + 1);
   limbs = (bits + FLINT_BITS - 1)/FLINT_BITS;

   TMP_START;

   rowm = (ulong *) TMP_ALLOC(n*sizeof(ulong)*limbs);
   flint_mpn_zero(rowm, n*limbs);

   for (i = 0, j = 0; i < n; i++, j += limbs)
      rowm[j] = A[m][i];
   
   for (i = 0; i < n; i++)
   {
      if (i != 0)
      {
         switch (limbs)
         {
            case 1:
               NMOD_RED(A[m][i], rowm[i], M->mod);
               break;
            case 2:
               NMOD2_RED2(A[m][i], rowm[2*i + 1], rowm[2*i], M->mod);
               break;
            case 3:
               NMOD_RED(rowm[3*i + 2], rowm[3*i + 2], M->mod);                                      \
               NMOD_RED3(A[m][i], rowm[3*i + 2], rowm[3*i + 1], rowm[3*i], M->mod);
               break;
         }
      }

      if (A[m][i] != 0)
      {
         r = P[i];
         if (r != -WORD(1))
         {
            h = n_negmod(A[m][i], M->mod.n);
            A[m][i] = 0;

            switch (limbs)
            {
            case 1:
               for (j = i + 1; j < L[r]; j++)
                  rowm[j] += A[r][j]*h;
               break;
            case 2:
               for (j = i + 1; j < L[r]; j++)
               {
                  umul_ppmm(hi, lo, A[r][j], h);
                  add_ssaaaa(rowm[2*j + 1], rowm[2*j], rowm[2*j + 1], rowm[2*j], hi, lo);
               }
               break;
            case 3:
               for (j = i + 1; j < L[r]; j++)
               {
                  umul_ppmm(hi, lo, A[r][j], h);
                  add_sssaaaaaa(rowm[3*j + 2], rowm[3*j + 1], rowm[3*j],
                                rowm[3*j + 2], rowm[3*j + 1], rowm[3*j], 0, hi, lo);
               }
               break;
            }
         } else
         {
            h = n_invmod(A[m][i], M->mod.n);
            A[m][i] = 1;

            switch (limbs)
            {
            case 1:
               for (j = i + 1; j < L[m]; j++)
               {
                  NMOD_RED(A[m][j], rowm[j], M->mod);
                  A[m][j] = n_mulmod2_preinv(A[m][j], h, M->mod.n, M->mod.ninv);
               }
               break;
            case 2:
               for (j = i + 1; j < L[m]; j++)
               {
                  NMOD2_RED2(A[m][j], rowm[2*j + 1], rowm[2*j], M->mod);
                  A[m][j] = n_mulmod2_preinv(A[m][j], h, M->mod.n, M->mod.ninv);
               }
               break;
            case 3:
               for (j = i + 1; j < L[m]; j++)
               {
                  NMOD_RED(rowm[3*j + 2], rowm[3*j + 2], M->mod);                                      \
                  NMOD_RED3(A[m][j], rowm[3*j + 2], rowm[3*j + 1], rowm[3*j], M->mod);
                  A[m][j] = n_mulmod2_preinv(A[m][j], h, M->mod.n, M->mod.ninv);
               }
               break;
            }

            P[i] = m;

            res = i;

            break;
         }
      }
   }

   TMP_END;

   return res;
}
