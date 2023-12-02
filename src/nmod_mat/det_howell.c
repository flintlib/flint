/*
    Copyright (C) 2015 Tommy Hofmann
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_mat.h"

/*
  We do no quite need the full Howell form. We just need to reduce to upper
  triangular form using row swaps and Euclidean row operations (i.e. that
  use the Euclidean xgcd to replace the pivot with its gcd with the entry
  in the row in question (mod n). Of course we need a special version of the
  Euclidean xgcd whose first cofactor is a unit mod n.
*/

/*
   Find s, t such that g = s*a - t*b is the gcd of a and b mod n and where
   s is a unit mod n. Assumes a and b are reduced mod n and no aliasing.
*/
static inline mp_limb_t
_nmod_xgcd_unit(mp_limb_t * s, mp_limb_t * t, mp_limb_t a, mp_limb_t b, nmod_t mod)
{
   mp_limb_t g, ag, bg;

   if (a >= b)
      g = n_xgcd(s, t, a, b);
   else /* b > a */
   {
      g = n_xgcd(t, s, b, a);
      *s = nmod_neg(*s, mod);
      *t = nmod_neg(*t, mod);
   }

   ag = a/g;
   bg = b/g;

   while (n_gcd(*s, mod.n) != 1)
   {
      *s = nmod_add(*s, bg, mod);
      *t = nmod_add(*t, ag, mod);
   }

   return g;
}

static inline int
_nmod_mat_pivot(nmod_mat_t A, slong start_row, slong col)
{
    slong j;
    mp_ptr u;

    if (nmod_mat_entry(A, start_row, col) != 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (nmod_mat_entry(A, j, col) != 0)
        {
            u = A->rows[j];
            A->rows[j] = A->rows[start_row];
            A->rows[start_row] = u;

            return -1;
        }
    }
    return 0;
}

/* test whether q*a = b mod N has a solution */
static int
_n_is_divisible(mp_ptr q, mp_limb_t b, mp_limb_t a, nmod_t N)
{
    mp_limb_t e, g;
    g = n_gcdinv(&e, a, N.n);

    if (( b % g ) == 0)
    {
        *q = nmod_mul(e, b/g, N);
        return 1;
    }

    return 0;
}

mp_limb_t _nmod_mat_det_howell(nmod_mat_t A)
{
    mp_limb_t s, t, t1, det = 1, unit = 1;
    slong m, n, row, col, i, k;
    nmod_t mod = A->mod;

    if (nmod_mat_is_empty(A))
        return mod.n != 1;

    n = A->r;
    m = A->c;

    row = col = 0;

    while (row < n && col < m)
    {
        int pivswap = _nmod_mat_pivot(A, row, col);

	if (pivswap == 0)
            return 0;

	if (pivswap == -1)
           det = nmod_neg(det, mod);

	for (i = row + 1; i < n; i++)
        {
            if (nmod_mat_entry(A, i, col) == 0)
                continue;

            if (_n_is_divisible(&s, nmod_mat_entry(A, i, col), nmod_mat_entry(A, row, col), mod))
            {
                 for (k = col; k < m; k++)
                 {
                     t = nmod_sub(nmod_mat_entry(A, i, k), nmod_mul(s, nmod_mat_entry(A, row, k), mod), mod);
                     nmod_mat_entry(A, i, k) = t;
                 }
            }
            else
            {
                _nmod_xgcd_unit(&s, &t, nmod_mat_entry(A, row, col), nmod_mat_entry(A, i, col), mod);

                /* now g = s*x - t*y mod n */
		unit = nmod_mul(unit, s, mod);

                for (k = col; k < m; k++)
                {

                    t1 = nmod_sub(nmod_mul(s, nmod_mat_entry(A, row, k), mod), nmod_mul(t, nmod_mat_entry(A, i, k), mod), mod);
                    nmod_mat_entry(A, row, k) = t1;
                }

		/* now it's divisible, restart this row */
                i--;
		continue;
            }
        }

        det = nmod_mul(det, nmod_mat_entry(A, row, col), mod);

        row++;
        col++;
    }

    unit = nmod_inv(unit, mod);

    return nmod_mul(det, unit, mod);
}

mp_limb_t
nmod_mat_det_howell(const nmod_mat_t A)
{
    nmod_mat_t tmp;
    mp_limb_t det;
    slong dim = A->r;

    if (dim != A->c)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_mat_det_howell). Non-square matrix.\n");
    }

    nmod_mat_init_set(tmp, A);
    det = _nmod_mat_det_howell(tmp);
    nmod_mat_clear(tmp);

    return det;
}
