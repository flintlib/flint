/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"

void fmpz_mat_hnf_mod_D(fmpz_mat_t H, const fmpz_mat_t A, const fmpz_t D)
{
    slong j, j2, i, k;
    fmpz_t R, R2;

    fmpz_init_set(R, D);
    fmpz_init(R2);
    fmpz_mat_set(H, A);
    for (j = 0, k = 0; j != A->c; j++, k++)
    {
        fmpz_t d, u, v;
        fmpz_init(u);
        fmpz_init(v);
        fmpz_init(d);
        fmpz_cdiv_q_ui(R2, R, 2);

        if (fmpz_is_zero(fmpz_mat_entry(H, k, j)))
            fmpz_set(fmpz_mat_entry(H, k, j), R);
        for (i = k + 1; i != A->r; i++)
        {
            /* reduce row i with row k mod R */
            fmpz_t r1d, r2d, b;
            if (fmpz_is_zero(fmpz_mat_entry(H, i, j)))
                continue;
            fmpz_init(r1d);
            fmpz_init(r2d);
            fmpz_init(b);
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, k, j), fmpz_mat_entry(H, i, j));
            fmpz_divexact(r2d, fmpz_mat_entry(H, i, j), d);
            fmpz_divexact(r1d, fmpz_mat_entry(H, k, j), d);
            for (j2 = 0; j2 < A->c; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, k, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, i, j2));
                fmpz_mul(fmpz_mat_entry(H, i, j2), r1d, fmpz_mat_entry(H, i, j2));
                fmpz_submul(fmpz_mat_entry(H, i, j2), r2d, fmpz_mat_entry(H, k, j2));
                fmpz_mod(fmpz_mat_entry(H, i, j2), fmpz_mat_entry(H, i, j2), R);
                if (fmpz_cmp(fmpz_mat_entry(H, i, j2), R2) > 0)
                    fmpz_sub(fmpz_mat_entry(H, i, j2), fmpz_mat_entry(H, i, j2), R);
                fmpz_set(fmpz_mat_entry(H, k, j2), b);
                fmpz_mod(fmpz_mat_entry(H, k, j2), fmpz_mat_entry(H, k, j2), R);
                if (fmpz_cmp(fmpz_mat_entry(H, k, j2), R2) > 0)
                    fmpz_sub(fmpz_mat_entry(H, k, j2), fmpz_mat_entry(H, k, j2), R);
            }
            fmpz_clear(b);
            fmpz_clear(r2d);
            fmpz_clear(r1d);
        }
        fmpz_xgcd(d, u, v, fmpz_mat_entry(H, k, j), R);
        for (j2 = 0; j2 < A->c; j2++)
        {
            fmpz_mul(fmpz_mat_entry(H, k, j2), u,
                fmpz_mat_entry(H, k, j2));
            fmpz_mod(fmpz_mat_entry(H, k, j2), fmpz_mat_entry(H, k, j2), R);
        }
        if (fmpz_is_zero(fmpz_mat_entry(H, k, j)))
            fmpz_set(fmpz_mat_entry(H, k, j), R);
        /* reduce higher entries of column j with row k */
        for (i = k - 1; i >= 0; i--)
        {
            fmpz_t q;
            fmpz_init(q);
            fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                    fmpz_mat_entry(H, k, j));
            for (j2 = 0; j2 < A->c; j2++)
            {
                fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                        fmpz_mat_entry(H, k, j2));
                fmpz_mod(fmpz_mat_entry(H, i, j2), fmpz_mat_entry(H, i, j2), D);
            }
            fmpz_clear(q);
        }
        fmpz_divexact(R, R, d);
        fmpz_clear(d);
        fmpz_clear(v);
        fmpz_clear(u);
    }
    fmpz_clear(R2);
    fmpz_clear(R);
}
