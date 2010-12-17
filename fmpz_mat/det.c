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


void _fmpz_mat_det_2x2(fmpz_t det, fmpz ** const x)
{
    /* Use a temporary in case someone gets the bright idea to
       save the determinant in one of the entries of
       the input matrix */
    /*
    fmpz_t t;
    fmpz_init(t);
    */

    fmpz_mul   (det, &x[0][0], &x[1][1]);
    fmpz_submul(det, &x[0][1], &x[1][0]);

    /*
    fmpz_set(det, t);
    fmpz_clear(t);
    */
}

void _fmpz_mat_det_3x3(fmpz_t det, fmpz ** const x)
{
    fmpz_t a, t;

    fmpz_init(a);
    fmpz_init(t);

    fmpz_mul   (a, &x[1][0], &x[2][1]);
    fmpz_submul(a, &x[1][1], &x[2][0]);
    fmpz_mul   (t, a, &x[0][2]);

    fmpz_mul   (a, &x[1][2], &x[2][0]);
    fmpz_submul(a, &x[1][0], &x[2][2]);
    fmpz_addmul(t, a, &x[0][1]);

    fmpz_mul   (a, &x[1][1], &x[2][2]);
    fmpz_submul(a, &x[1][2], &x[2][1]);
    fmpz_addmul(t, a, &x[0][0]);

    fmpz_set(det, t);

    fmpz_clear(a);
    fmpz_clear(t);
}

void _fmpz_mat_det_4x4(fmpz_t det, fmpz ** const x)
{
    fmpz_t a, b, t;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(t);

    fmpz_mul   (a, &x[0][3], &x[1][2]);
    fmpz_submul(a, &x[0][2], &x[1][3]);
    fmpz_mul   (b, &x[2][1], &x[3][0]);
    fmpz_submul(b, &x[2][0], &x[3][1]);
    fmpz_mul(t, a, b);

    fmpz_mul   (a, &x[0][1], &x[1][3]);
    fmpz_submul(a, &x[0][3], &x[1][1]);
    fmpz_mul   (b, &x[2][2], &x[3][0]);
    fmpz_submul(b, &x[2][0], &x[3][2]);
    fmpz_addmul(t, a, b);

    fmpz_mul   (a, &x[0][2], &x[1][1]);
    fmpz_submul(a, &x[0][1], &x[1][2]);
    fmpz_mul   (b, &x[2][3], &x[3][0]);
    fmpz_submul(b, &x[2][0], &x[3][3]);
    fmpz_addmul(t, a, b);

    fmpz_mul   (a, &x[0][3], &x[1][0]);
    fmpz_submul(a, &x[0][0], &x[1][3]);
    fmpz_mul   (b, &x[2][2], &x[3][1]);
    fmpz_submul(b, &x[2][1], &x[3][2]);
    fmpz_addmul(t, a, b);

    fmpz_mul   (a, &x[0][0], &x[1][2]);
    fmpz_submul(a, &x[0][2], &x[1][0]);
    fmpz_mul   (b, &x[2][3], &x[3][1]);
    fmpz_submul(b, &x[2][1], &x[3][3]);
    fmpz_addmul(t, a, b);

    fmpz_mul   (a, &x[0][1], &x[1][0]);
    fmpz_submul(a, &x[0][0], &x[1][1]);
    fmpz_mul   (b, &x[2][3], &x[3][2]);
    fmpz_submul(b, &x[2][2], &x[3][3]);
    fmpz_addmul(t, a, b);

    fmpz_set(det, t);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(t);
}

void
_fmpz_mat_det_rowreduce(fmpz_t det, const fmpz_mat_t A)
{
    fmpz_mat_t tmp;
    long m = A->r;
    long rank;

    fmpz_mat_init_set(tmp, A);
    rank = _fmpz_mat_rowreduce(tmp, ROWREDUCE_FAST_ABORT);

    if (rank < 0)
    {
        rank = -rank;
        if (rank < m)
            fmpz_zero(det);
        else
            fmpz_neg(det, &tmp->rows[rank-1][rank-1]);
    }
    else
    {
        if (rank < m)
            fmpz_zero(det);
        else
            fmpz_set(det, &tmp->rows[rank-1][rank-1]);
    }
    fmpz_mat_clear(tmp);
}

void
fmpz_mat_det(fmpz_t det, const fmpz_mat_t A)
{
    long dim = A->r;

    if (dim != A->c)
    {
        printf("fmpz_mat_det: nonsquare matrix");
        abort();
    }

    switch (dim)
    {
        case 0:  fmpz_set_ui(det, 1UL);           break;
        case 1:  fmpz_set(det, A->entries);       break;
        case 2:  _fmpz_mat_det_2x2(det, A->rows); break;
        case 3:  _fmpz_mat_det_3x3(det, A->rows); break;
        case 4:  _fmpz_mat_det_4x4(det, A->rows); break;
        default: _fmpz_mat_det_rowreduce(det, A);
    }
}
