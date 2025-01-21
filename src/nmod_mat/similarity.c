/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_mat.h"

void nmod_mat_similarity(nmod_mat_t M, slong r, ulong d)
{
   slong n = M->r, i, j;

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
         NMOD_ADDMUL(nmod_mat_entry(M, i, j), nmod_mat_entry(M, i, r), d, M->mod);

      for (j = r + 1; j < n; j++)
         NMOD_ADDMUL(nmod_mat_entry(M, i, j), nmod_mat_entry(M, i, r), d, M->mod);
   }

   d = n_negmod(d, M->mod.n);

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
         NMOD_ADDMUL(nmod_mat_entry(M, r, i), nmod_mat_entry(M, j, i), d, M->mod);

      for (j = r + 1; j < n; j++)
         NMOD_ADDMUL(nmod_mat_entry(M, r, i), nmod_mat_entry(M, j, i), d, M->mod);
   }
}
