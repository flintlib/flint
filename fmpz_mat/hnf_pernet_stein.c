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
#include "fmpq_mat.h"

void add_column(fmpz_mat_t H, const fmpz_mat_t B, const fmpz_mat_t H1)
{
    slong i, j;
    fmpz_t M;
    fmpz_mat_t Bu, col, x;
    fmpq_mat_t x_q, col_q;

    fmpz_init(M);
    fmpz_mat_init(Bu, B->r, B->c - 1);
    fmpz_mat_init(col, B->r, 1);
    fmpz_mat_init(x, B->r, 1);
    fmpq_mat_init(x_q, B->r, 1);
    fmpq_mat_init(col_q, B->r, 1);

    for (i = 0; i < B->r; i++)
        fmpz_set(fmpz_mat_entry(col, i, 0), fmpz_mat_entry(B, i, B->c - 1));
    /* TODO: use random vector u */
    for (i = 0; i < B->r; i++)
        for (j = 0; j < B->c - 1; j++)
            fmpz_set(fmpz_mat_entry(Bu, i, j), fmpz_mat_entry(B, i, j));
    fmpz_mat_solve_dixon(x, M, Bu, col);
    fmpq_mat_set_fmpz_mat_mod_fmpz(x_q, x, M);
    fmpz_mat_mul(col, H1, x);
    fmpq_mat_set_fmpz_mat_mod_fmpz(col_q, col, M);
    fmpq_mat_get_fmpz_mat(col, col_q);
    for (i = 0; i < H1->r; i++)
    {
        for (j = 0; j < H1->c; j++)
            fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H1, i, j));
        fmpz_set(fmpz_mat_entry(H, i, H1->c), fmpz_mat_entry(col, i, 0));
    }

    fmpq_mat_clear(x_q);
    fmpq_mat_clear(col_q);
    fmpz_mat_clear(x);
    fmpz_mat_clear(col);
    fmpz_mat_clear(Bu);
    fmpz_clear(M);
}

void add_row(fmpz_mat_t H, const fmpz_mat_t A, const fmpz * row)
{
    slong i, i2, j, j2, num_pivots, new_row;
    slong * pivots;
    fmpz_t b, d, u, v, r1d, r2d, q;

    fmpz_init(b);
    fmpz_init(d);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(r1d);
    fmpz_init(r2d);
    fmpz_init(q);
    num_pivots = 0;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(A, i, j));
    for (j = 0; j < A->c; j++)
        fmpz_set(fmpz_mat_entry(H, i, j), &row[j]);

    /* find the pivots of A */
    for (i = j = 0; i < A->r; i++, j++)
    {
        for (; j < A->c && fmpz_is_zero(fmpz_mat_entry(A, i, j)); j++);
        if (j == A->c)
            break;
        num_pivots = i + 1;
    }
    pivots = flint_calloc(num_pivots, sizeof(slong));
    for (i = j = 0; i < A->r; i++, j++)
    {
        for (; j < A->c && fmpz_is_zero(fmpz_mat_entry(A, i, j)); j++);
        if (j == A->c    )
            break;
        pivots[i] = j;
    }

    /* reduce row to be added with existing */
    for (i = 0; i < num_pivots; i++)
    {
        j = pivots[i];
        if (fmpz_is_zero(fmpz_mat_entry(H, A->r, j)))
            continue;
        fmpz_xgcd(d, u, v, fmpz_mat_entry(H, i, j), fmpz_mat_entry(H, A->r, j));
        fmpz_divexact(r1d, fmpz_mat_entry(H, i, j), d);
        fmpz_divexact(r2d, fmpz_mat_entry(H, A->r, j), d);
        for (j2 = j; j2 < A->c; j2++)
        {
            fmpz_mul(b, u, fmpz_mat_entry(H, i, j2));
            fmpz_addmul(b, v, fmpz_mat_entry(H, A->r, j2));
            fmpz_mul(fmpz_mat_entry(H, A->r, j2), r1d, fmpz_mat_entry(H, A->r, j2));
            fmpz_submul(fmpz_mat_entry(H, A->r, j2), r2d, fmpz_mat_entry(H, i, j2));
            fmpz_set(fmpz_mat_entry(H, i, j2), b);
        }
    }
    /* find first non-zero entry of the added row */
    for (j = 0; j < A->c && fmpz_is_zero(fmpz_mat_entry(H, A->r, j)); j++);
    new_row = A->r;
    if (j != A->c) /* last row non-zero, move to correct position */
    {
        if (fmpz_sgn(fmpz_mat_entry(H, A->r, j)) < 0)
        {
            for (j2 = j; j2 < A->c; j2++)
            {
                fmpz_neg(fmpz_mat_entry(H, A->r, j2),
                        fmpz_mat_entry(H, A->r, j2));
            }
        }
        do
        {
            if (new_row < A->r)
                fmpz_mat_swap_rows(H, NULL, new_row, new_row + 1);
            new_row--;
            for (j2 = 0; j2 < A->c && fmpz_is_zero(fmpz_mat_entry(H, new_row, j2)); j2++);
        }
        while (j2 > j);
    }

    for (i = j = num_pivots = 0; i < H->r; i++, j++)
    {
        for (; j < H->c && fmpz_is_zero(fmpz_mat_entry(H, i, j)); j++);
        if (j == H->c)
            break;
        num_pivots = i + 1;
    }
    pivots = flint_realloc(pivots, num_pivots * sizeof(slong));
    for (i = new_row, j = 0; i < H->r; i++, j++)
    {
        for (; j < H->c && fmpz_is_zero(fmpz_mat_entry(H, i, j)); j++);
        if (j == H->c)
            break;
        pivots[i] = j;
    }

    for (i = 0; i < num_pivots; i++)
    {
        for (i2 = 0; i2 < i; i2++)
        {
            fmpz_fdiv_q(q, fmpz_mat_entry(H, i2, pivots[i]),
                    fmpz_mat_entry(H, i, pivots[i]));
            for (j2 = pivots[i]; j2 < A->c; j2++)
            {
                fmpz_submul(fmpz_mat_entry(H, i2, j2), q,
                        fmpz_mat_entry(H, i, j2));
            }
        }
    }

    fmpz_clear(q);
    fmpz_clear(r2d);
    fmpz_clear(r1d);
    fmpz_clear(v);
    fmpz_clear(u);
    fmpz_clear(d);
    fmpz_clear(b);
    flint_free(pivots);
}

void fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong i, j;
    fmpz_t d1, d2, g, s, t;
    fmpz_mat_t B, C, H1, H2, H3, H4;

    fmpz_init(d1);
    fmpz_init(d2);
    fmpz_init(g);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_mat_init(B, A->r - 1, A->c);
    fmpz_mat_init(C, A->r - 1, A->c - 1);
    fmpz_mat_init(H1, A->r - 1, A->c - 1);
    fmpz_mat_init(H2, A->r - 1, A->c);
    fmpz_mat_init(H3, A->r, A->c);
    fmpz_mat_init(H4, A->r + 1, A->c);
    for (i = 0; i < A->r - 1; i++)
    {
        for (j = 0; j < A->c - 1; j++)
        {
            fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j));
            fmpz_set(fmpz_mat_entry(C, i, j), fmpz_mat_entry(A, i, j));
        }
        fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j));
    }
    fmpz_mat_det(d1, C);
    for (j = 0; j < A->c - 1; j++)
        fmpz_set(fmpz_mat_entry(C, A->r - 2, j), fmpz_mat_entry(A, A->r - 1, j));
    fmpz_mat_det(d2, C);
    fmpz_xgcd(g, s, t, d1, d2);
    for (j = 0; j < A->c - 1; j++)
    {
        fmpz_mul(fmpz_mat_entry(C, A->r - 2, j), s, fmpz_mat_entry(A, A->r - 2, j));
        fmpz_addmul(fmpz_mat_entry(C, A->r - 2, j), t, fmpz_mat_entry(A, A->r - 1, j));
        fmpz_set(fmpz_mat_entry(B, A->r - 2, j), fmpz_mat_entry(C, A->r - 2, j));
    }
    fmpz_mul(fmpz_mat_entry(B, A->r - 2, A->c - 1), s, fmpz_mat_entry(A, A->r - 2, A->c - 1));
    fmpz_addmul(fmpz_mat_entry(B, A->r - 2, A->c - 1), t, fmpz_mat_entry(A, A->r - 1, A->c - 1));
    if (!fmpz_is_zero(g)) /* chosen matrix non-invertible */
    {
        /* if g is too big, recurse */
        fmpz_abs(g, g);
        if (COEFF_IS_MPZ(*g) && 0)
        {
            fmpz_print(g);
            flint_printf("\n\n");
        }
        if (COEFF_IS_MPZ(*g) && C->r > 3)
        {
            fmpz_mat_hnf_pernet_stein(H1, C);
        }
        else /* use modulo determinant algorithm to compute HNF of C */
        {
            fmpz_mat_hnf_mod_D(H1, C, g);
        }
        add_column(H2, B, H1);
        add_row(H3, H2, A->rows[A->r - 2]);
        add_row(H4, H3, A->rows[A->r - 1]);
        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H4, i, j));
    }
    else
    {
        fmpz_mat_hnf_xgcd(H, A);
    }
    fmpz_clear(d2);
    fmpz_clear(d1);
    fmpz_clear(t);
    fmpz_clear(s);
    fmpz_clear(g);
    fmpz_mat_clear(H1);
    fmpz_mat_clear(H2);
    fmpz_mat_clear(H3);
    fmpz_mat_clear(H4);
    fmpz_mat_clear(C);
    fmpz_mat_clear(B);
}