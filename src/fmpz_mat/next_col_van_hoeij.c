/*
    Copyright (C) 2011, 2016 William Hart
    Copyright (C) 2011 Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <gmp.h>
#include "longlong.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

/* todo: proper inplace realloc */
static void _fmpz_mat_resize_van_hoeij(fmpz_mat_t M, slong r, slong c)
{
    slong i;
    fmpz_mat_t T;

    FLINT_ASSERT(r >= M->r);
    FLINT_ASSERT(c >= M->c);

    fmpz_mat_init(T, r, c);
    for (i = 0; i < M->r; i++)
        _fmpz_vec_swap(fmpz_mat_row(T, i + r - M->r), fmpz_mat_row(M, i), M->c);
    fmpz_mat_swap(M, T);
    fmpz_mat_clear(T);
}

/*
    Goal here is to take a matrix M, get U, multiply U by col look at max bits
    of U*col and P subtract exp and decide if it's worth calling LLL.
    If is not return 0. U_exp is the assumed power of the scalar multiple of U
    (so 2^U_exp is the scalar weight of U).
    This should match the virtual precision.
    If it is then return weight of last column (check that it should not
    be zero... ??) and  augment M with the new column and new row
    This new column should be truncated to the correct amount before
    re-multiplying by U first make sure there are enough bits to even bother.
*/
int fmpz_mat_next_col_van_hoeij(fmpz_mat_t M, fmpz_t P,
                                        fmpz_mat_t col, slong exp, slong U_exp)
{
   slong j, k, r = col->r;
   slong bit_r = FLINT_MAX(r, 20);
   slong s = M->r;
   fmpz_mat_t U, x, y;
   fmpz_t P_trunc;
   slong mbts_full;

   k = fmpz_bits(P) - bit_r - bit_r/2;

   /* primary check: is LLL potentially justified at this precision? */
   if (k < exp + (slong) FLINT_BIT_COUNT(r + 1))
      return 0;

   fmpz_init(P_trunc);
   fmpz_mat_init(x, r, 1);
   fmpz_mat_init(y, s, 1);

   /* find U, the combinatorial part of M */
   fmpz_mat_window_init(U, M, 0, 0, s, r);

   /*
      Secondary check (from FLINT 1.6, _F_mpz_mat_next_col):
      compute y_full = U * col / 2^U_exp mod P at full CLD precision (no
      k-based truncation).  If the sup-norm of y_full is too small relative
      to the CLD bound (exp), the column carries insufficient information to
      improve the lattice even though the primary ISD check passed.  This
      filters out columns whose actual data is negligibly small -- typically
      because the CLD coefficient is genuinely close to a multiple of P at
      this precision -- and avoids a full LLL call that would remove zero rows.

      Threshold 0.973*bit_r matches FLINT 1.6.
   */
   fmpz_mat_t y_full;
   fmpz_mat_init(y_full, s, 1);
   fmpz_mat_mul(y_full, U, col);
   fmpz_mat_scalar_tdiv_q_2exp(y_full, y_full, U_exp);
   fmpz_mat_scalar_smod(y_full, y_full, P);
   mbts_full = FLINT_ABS(fmpz_mat_max_bits(y_full));

   if (mbts_full < (slong)(0.973*(double)bit_r - 0.1) + exp)
   {
      fmpz_mat_clear(x);
      fmpz_mat_clear(y);
      fmpz_mat_clear(y_full);
      fmpz_clear(P_trunc);
      fmpz_mat_window_clear(U);
      return 0;
   }

   fmpz_mat_clear(y_full);

   /* find scale factor 2^k */
   k -= U_exp; /* we want this many bits beyond the radix point */

   if (k >= 0)
   {
      fmpz_mat_scalar_tdiv_q_2exp(x, col, k);
      fmpz_tdiv_q_2exp(P_trunc, P, k);
   } else
   {
      fmpz_mat_scalar_mul_2exp(x, col, -k);
      fmpz_mul_2exp(P_trunc, P, -k);
   }

   /* multiply column by U */
   fmpz_mat_mul(y, U, x);

   /* everything in U was already scaled by U_exp, so divide out scaling */
   fmpz_mat_scalar_tdiv_q_2exp(y, y, U_exp);
   fmpz_mat_scalar_smod(y, y, P_trunc);

   /* resize M */
   _fmpz_mat_resize_van_hoeij(M, s + 1, M->c + 1);

   /* insert new column and row data */
   fmpz_set(fmpz_mat_entry(M, 0, M->c - 1), P_trunc);

   for (j = 1; j < M->r; j++)
      fmpz_set(fmpz_mat_entry(M, j, M->c - 1), fmpz_mat_entry(y, j - 1, 0));

   fmpz_mat_clear(x);
   fmpz_mat_clear(y);
   fmpz_clear(P_trunc);

   fmpz_mat_window_clear(U);

   return 1;
}
