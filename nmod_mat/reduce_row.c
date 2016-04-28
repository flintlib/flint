/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

slong nmod_mat_reduce_row(nmod_mat_t M, slong * P, slong * L, slong m)
{
   slong n = M->c, i, j, r;
   ulong ** A = M->rows;
   ulong h;

   for (i = 0; i < n; i++)
   {
      if (A[m][i] != 0)
      {
         r = P[i];
         if (r != -WORD(1))
         {
            h = n_negmod(A[m][i], M->mod.n);
            A[m][i] = 0;

            for (j = i + 1; j < L[r]; j++)
               NMOD_ADDMUL(A[m][j], A[r][j], h, M->mod);
         } else
         {
            h = n_invmod(A[m][i], M->mod.n);
            A[m][i] = 1;

            for (j = i + 1; j < L[m]; j++)
               A[m][j] = n_mulmod2_preinv(A[m][j], h, M->mod.n, M->mod.ninv);

            P[i] = m;

            return i;
         }
      }
   }
   
   return -WORD(1);
}
