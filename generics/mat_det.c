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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"
#include "perm.h"

long _elem_mat_fflu(elem_ptr * mat, elem_ptr den, long * perm, long m, long n, int rank_check, const ring_t ring);


/* Computes the determinant using Gaussian elimination,
   destroying the input matrix */

void
_elem_mat_det(elem_ptr det, elem_ptr * mat, long n, const ring_t ring)
{
    long * perm = _perm_init(n);

    _elem_mat_fflu(mat, det, perm, n, n, 1, ring);

    if (_perm_parity(perm, n))
        elem_neg(det, det, ring);

    _perm_clear(perm);
}

void
elem_mat_det(elem_ptr det, const elem_mat_t mat, const ring_t ring)
{
    long n = mat->r;
    const ring_struct * det_ring = RING_PARENT(ring);

    if (n != mat->c)
    {
        printf("mat_det: matrix must be square\n");
        abort();
    }

    if (n == 0)
    {
        elem_one(det, RING_PARENT(ring));
    }
    else if (n == 1)
    {
        elem_set(det, MAT_SRCINDEX(mat->rows, 0, 0, det_ring), det_ring);
    }
    else
    {
        elem_mat_t tmp;
        elem_mat_init(tmp, mat->r, mat->c, ring);
        elem_mat_set(tmp, mat, ring);
        _elem_mat_det(det, tmp->rows, n, det_ring);
        elem_mat_clear(tmp, ring);
    }
}

