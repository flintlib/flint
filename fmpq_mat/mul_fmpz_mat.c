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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"


void
fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, const fmpq_mat_t A, const fmpz_mat_t B)
{
    long i, j;

    fmpz_mat_t Aclear;
    fmpz_mat_t Cclear;

    fmpz * Aden;

    fmpz_mat_init(Aclear, A->r, A->c);
    fmpz_mat_init(Cclear, A->r, B->c);

    Aden = _fmpz_vec_init(A->r);

    fmpq_mat_get_fmpz_mat_rowwise(Aclear, Aden, A);
    fmpz_mat_mul(Cclear, Aclear, B);

    for (i = 0; i < C->r; i++)
    {
        for (j = 0; j < C->c; j++)
        {
            fmpz_set(fmpq_mat_entry_num(C, i, j), fmpz_mat_entry(Cclear, i, j));
            fmpz_set(fmpq_mat_entry_den(C, i, j), Aden + i);
            fmpq_canonicalise(fmpq_mat_entry(C, i, j));
        }
    }

    fmpz_mat_clear(Aclear);
    fmpz_mat_clear(Cclear);

    _fmpz_vec_clear(Aden, A->r);
}
