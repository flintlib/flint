/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "perm.h"
#include "longlong.h"

#define E(j,k) fmpz_mat_entry(B,j,k)

slong
fmpz_mat_fflu(fmpz_mat_t B, fmpz_t den, slong * perm,
                            const fmpz_mat_t A, int rank_check)
{
    slong m, n, j, k, rank, r, pivot_row, pivot_col, norm = 0;
    ulong p1h, p1l, p2h, p2l, uden = 0, dinv = 0, quo;
    ulong FLINT_SET_BUT_UNUSED(rem);
    slong mbits = fmpz_mat_max_bits(A);
    int small = FLINT_ABS(mbits) <= FLINT_BITS - 2;
    int dsgn = 0, sgn, den1 = 0, work_to_do;

    if (fmpz_mat_is_empty(A))
    {
        fmpz_one(den);
        return 0;
    }

    fmpz_mat_set(B, A);
    m = B->r;
    n = B->c;
    rank = pivot_row = pivot_col = 0;

    while (pivot_row < m && pivot_col < n)
    {
        r = fmpz_mat_find_pivot_any(B, pivot_row, m, pivot_col);

        if (r == -1)
        {
            if (rank_check)
            {
                fmpz_zero(den);
                rank = 0;
                break;
            }
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
            fmpz_mat_swap_rows(B, perm, pivot_row, r);

        rank++;

        if (small)
        {
            for (j = pivot_row + 1; j < m; j++)
            {
                work_to_do = !den1 || !fmpz_is_zero(E(j, pivot_col)) ||
                   !fmpz_is_one(E(pivot_row, pivot_col));

                if (work_to_do)
                {
                    for (k = pivot_col + 1; k < n; k++)
                    {
                        smul_ppmm(p1h, p1l, *E(j, k), *E(pivot_row, pivot_col));
                        smul_ppmm(p2h, p2l, *E(j, pivot_col), *E(pivot_row, k));
                        sub_ddmmss(p1h, p1l, p1h, p1l, p2h, p2l);

                        sgn = 0 > (slong) p1h;

                        if (sgn) /* take absolute value */
                           sub_ddmmss(p1h, p1l, UWORD(0), UWORD(0), p1h, p1l);

                        if (pivot_row > 0 && !den1)
                        {
                            if (p1h >= uden)
                            {
                                fmpz_set_uiui(E(j, k), p1h, p1l);

                                if (sgn)
                                    fmpz_neg(E(j, k), E(j, k));

                                fmpz_divexact(E(j, k), E(j, k), den);

                                small = 0;
                            } else
                            {
                                udiv_qrnnd_preinv(quo, rem,
                                  (p1h << norm) +
                                  r_shift(p1l, (FLINT_BITS - norm)),
                                      p1l << norm, uden << norm, dinv);

                                if (sgn ^ dsgn)
                                    fmpz_neg_ui(E(j, k), quo);
                                else
                                    fmpz_set_ui(E(j, k), quo);

                                if (quo > COEFF_MAX)
                                    small = 0;
                            }
                        } else
                        {
                            if (p1h > 0)
                            {
                                fmpz_set_uiui(E(j, k), p1h, p1l);

                                small = 0;
                            } else
                            {
                                fmpz_set_ui(E(j, k), p1l);

                                if (p1l > COEFF_MAX)
                                    small = 0;
                            }

                            if (sgn)
                                fmpz_neg(E(j, k), E(j, k));
                        }
                    }
                }
            }
        } else
        {
            for (j = pivot_row + 1; j < m; j++)
            {
                for (k = pivot_col + 1; k < n; k++)
                {
                    fmpz_mul(E(j, k), E(j, k), E(pivot_row, pivot_col));
                    fmpz_submul(E(j, k), E(j, pivot_col), E(pivot_row, k));

                    if (pivot_row > 0 && !den1)
                        fmpz_divexact(E(j, k), E(j, k), den);
                }
            }
        }

        fmpz_set(den, E(pivot_row, pivot_col));

        den1 = fmpz_is_one(den);

        if (small)
        {
            uden = FLINT_ABS((slong)(*den));
            dsgn = 0 > (slong)(*den);
            count_leading_zeros(norm, uden);
            invert_limb(dinv, uden << norm);

            if (fmpz_sizeinbase(den, 2) > FLINT_BITS - 2)
                small = 0;
        }

        pivot_row++;
        pivot_col++;
    }

    return rank;
}
