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

len_t
fmpq_mat_rref_fraction_free(fmpq_mat_t B, const fmpq_mat_t A)
{
    fmpz_mat_t Aclear;
    fmpz_t den;
    len_t rank;

    if (fmpq_mat_is_empty(A))
        return 0;

    fmpz_mat_init(Aclear, A->r, A->c);
    fmpq_mat_get_fmpz_mat_rowwise(Aclear, NULL, A);
    fmpz_init(den);

    rank = fmpz_mat_rref(Aclear, den, Aclear);

    if (rank == 0)
        fmpq_mat_zero(B);
    else
        fmpq_mat_set_fmpz_mat_div_fmpz(B, Aclear, den);

    fmpz_mat_clear(Aclear);
    fmpz_clear(den);

    return rank;
}
