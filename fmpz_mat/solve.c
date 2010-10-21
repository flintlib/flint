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

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"


void _fmpz_mat_solve_2x2(fmpz * x, fmpz_t d, fmpz ** const a, const fmpz * b)
{
    fmpz_mul   (&x[0], &a[1][1], &b[0]);
    fmpz_submul(&x[0], &a[0][1], &b[1]);
    fmpz_mul   (&x[1], &a[0][0], &b[1]);
    fmpz_submul(&x[1], &a[1][0], &b[0]);
    _fmpz_mat_det_2x2(d, a);
}


void _fmpz_mat_solve_3x3(fmpz * x, fmpz_t d, fmpz ** const a, const fmpz * b)
{
    fmpz_t t12, t13, t14, t15, t16, t17;

    fmpz_init(t12);
    fmpz_init(t13);
    fmpz_init(t14);
    fmpz_init(t15);
    fmpz_init(t16);
    fmpz_init(t17);

    fmpz_mul(t17, &a[1][0], &a[2][1]); fmpz_submul(t17, &a[1][1], &a[2][0]);
    fmpz_mul(t16, &a[1][2], &a[2][0]); fmpz_submul(t16, &a[1][0], &a[2][2]);
    fmpz_mul(t15, &a[1][1], &a[2][2]); fmpz_submul(t15, &a[1][2], &a[2][1]);
    fmpz_mul(t14, &a[2][0], &b[1]); fmpz_submul(t14, &a[1][0], &b[2]);
    fmpz_mul(t13, &a[2][1], &b[1]); fmpz_submul(t13, &a[1][1], &b[2]);
    fmpz_mul(t12, &a[2][2], &b[1]); fmpz_submul(t12, &a[1][2], &b[2]);

    fmpz_mul   (&x[0], t15, &b[0]);
    fmpz_addmul(&x[0], t13, &a[0][2]);
    fmpz_submul(&x[0], t12, &a[0][1]);

    fmpz_mul   (&x[1], t16, &b[0]);
    fmpz_addmul(&x[1], t12, &a[0][0]);
    fmpz_submul(&x[1], t14, &a[0][2]);

    fmpz_mul   (&x[2], t17, &b[0]);
    fmpz_addmul(&x[2], t14, &a[0][1]);
    fmpz_submul(&x[2], t13, &a[0][0]);

    fmpz_mul   (d, t15, &a[0][0]);
    fmpz_addmul(d, t16, &a[0][1]);
    fmpz_addmul(d, t17, &a[0][2]);

    fmpz_clear(t12);
    fmpz_clear(t13);
    fmpz_clear(t14);
    fmpz_clear(t15);
    fmpz_clear(t16);
    fmpz_clear(t17);
}

void
_fmpz_mat_solve_fflu(fmpz * x, fmpz_t den, const fmpz_mat_t A, const fmpz * b)
{
    long i, dim, rank;
    fmpz_mat_t T;
    fmpz * tmp;

    dim = A->r;

    /* Compute LU decomposition in a temporary matrix */
    fmpz_mat_init(T, dim, dim);
    _fmpz_vec_copy(T->entries, A->entries, dim * dim);
    rank = _fmpz_mat_rowreduce(T->rows, dim, dim, 0);

    if (FLINT_ABS(rank) == dim)
    {
        fmpz_set(den, &T->rows[dim-1][dim-1]);

        tmp = _fmpz_vec_init(dim);
        for (i = 0; i < dim; i++)
        {
            /* Insert in same order as the pivots */
            fmpz_set(&tmp[i], &b[(T->rows[i] - T->entries) / dim]);
        }

        _fmpz_mat_solve_fflu_precomp(tmp, T->rows, dim);

        for (i = 0; i < dim; i++)
            fmpz_set(&x[i], &tmp[i]);

        _fmpz_vec_clear(tmp, dim);
    }
    else
    {
        fmpz_zero(den);
    }

    fmpz_mat_clear(T);
}

void
fmpz_mat_solve(fmpz * x, fmpz_t den, const fmpz_mat_t A, const fmpz * b)
{
    long dim = A->r;

    FMPZ_MAT_ASSERT(dim == A->c, "fmpz_mat_solve: matrix must be square");

    switch (dim)
    {
        case 0:
            break;
        case 1:
            fmpz_set(den, A->entries);
            fmpz_set(x, b);
            break;
        case 2:
            _fmpz_mat_solve_2x2(x, den, A->rows, b);
            break;
        case 3:
            _fmpz_mat_solve_3x3(x, den, A->rows, b);
            break;
        default:
            _fmpz_mat_solve_fflu(x, den, A, b);
    }
}
