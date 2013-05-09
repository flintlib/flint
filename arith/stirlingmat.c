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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "arith.h"


void
_arith_stirling_next_row(fmpz * new, fmpz * prev, len_t n, len_t klen, int kind)
{
    len_t k;
    fmpz_t t, u;

    if (n == 0)
    {
        fmpz_one(new);
        return;
    }

    if (klen <= 0)
        return;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_zero(new);

    if (klen > n)
        fmpz_one(new + n);

    for (k = 1; k < FLINT_MIN(n, klen); k++)
    {
        fmpz_set(u, prev + k);
        fmpz_set(new + k, t);
        switch (kind)
        {
        case 0:
            fmpz_addmul_ui(new + k, u, n - 1UL);
            break;
        case 1:
            fmpz_submul_ui(new + k, u, n - 1UL);
            break;
        case 2:
            fmpz_addmul_ui(new + k, u, k);
            break;
        }
        fmpz_set(t, u);
    }
    fmpz_clear(t);
    fmpz_clear(u);
}

static void
_arith_stirling_mat(fmpz ** rows, len_t r, len_t c, int kind)
{
    len_t i, j;

    if (r == 0 || c == 0)
        return;

    fmpz_one(rows[0]);
    for (i = 1; i < c; i++)
        fmpz_zero(rows[0] + i);

    for (i = 1; i < r; i++)
    {
        _arith_stirling_next_row(rows[i], rows[i-1], i, FLINT_MIN(c,i+1), kind);

        for (j = i + 1; j < c; j++)
            fmpz_zero(rows[i] + j);
    }
}

void arith_stirling_number_1u_vec_next(fmpz * row, fmpz * prev, len_t n, len_t klen)
{
    _arith_stirling_next_row(row, prev, n, klen, 0);
}

void arith_stirling_number_1_vec_next(fmpz * row, fmpz * prev, len_t n, len_t klen)
{
    _arith_stirling_next_row(row, prev, n, klen, 1);
}

void arith_stirling_number_2_vec_next(fmpz * row, fmpz * prev, len_t n, len_t klen)
{
    _arith_stirling_next_row(row, prev, n, klen, 2);
}

void
arith_stirling_matrix_1u(fmpz_mat_t mat)
{
    _arith_stirling_mat(mat->rows, mat->r, mat->c, 0);
}

void
arith_stirling_matrix_1(fmpz_mat_t mat)
{
    _arith_stirling_mat(mat->rows, mat->r, mat->c, 1);
}

void
arith_stirling_matrix_2(fmpz_mat_t mat)
{
    _arith_stirling_mat(mat->rows, mat->r, mat->c, 2);
}
