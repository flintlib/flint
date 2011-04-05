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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A)
{
    long m, n, i, j, rank;

    fmpz_mat_t Aclear;
    fmpz * den;

    m = A->r;
    n = A->c;

    if (m <= 0 || n <= 0)
        return;

    fmpz_mat_init(Aclear, m, n);

    den = _fmpz_vec_init(m);

    fmpq_mat_get_fmpz_mat_rowwise(Aclear, den, A);
    rank = _fmpz_mat_rowreduce(NULL, Aclear, ROWREDUCE_FULL);

    if (rank == 0)
        fmpz_set_ui(den, 1UL); /* XXX: could just zero matrix */
    else
        fmpz_set(den, fmpz_mat_entry(Aclear, 0, 0));

    for (i = 0; i < B->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            fmpz_set(fmpq_mat_entry_num(B, i, j), fmpz_mat_entry(Aclear, i, j));
            fmpz_set(fmpq_mat_entry_den(B, i, j), den);
            fmpq_canonicalise(fmpq_mat_entry(B, i, j));
        }
    }

    fmpz_mat_clear(Aclear);

    _fmpz_vec_clear(den, m);
}
