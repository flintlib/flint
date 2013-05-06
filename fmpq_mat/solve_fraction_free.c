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

int
fmpq_mat_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A,
                                                    const fmpq_mat_t B)
{
    fmpz_mat_t Anum;
    fmpz_mat_t Bnum;
    fmpz_mat_t Xnum;
    fmpz_t den;
    int success;

    fmpz_mat_init(Anum, A->r, A->c);
    fmpz_mat_init(Bnum, B->r, B->c);
    fmpz_mat_init(Xnum, B->r, B->c);
    fmpz_init(den);

    fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
    success = fmpz_mat_solve(Xnum, den, Anum, Bnum);

    if (success)
        fmpq_mat_set_fmpz_mat_div_fmpz(X, Xnum, den);

    fmpz_mat_clear(Anum);
    fmpz_mat_clear(Bnum);
    fmpz_mat_clear(Xnum);
    fmpz_clear(den);

    return success;
}
