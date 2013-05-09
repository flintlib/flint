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

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "nmod_vec.h"


void
fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A)
{
    fmpz_t p, s, t;
    len_t i, j;

    fmpz_init(p);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_one(p);

    for (i = 0; i < A->r; i++)
    {
        fmpz_zero(s);

        for (j = 0; j < A->c; j++)
            fmpz_addmul(s, A->rows[i] + j, A->rows[i] + j);

        fmpz_sqrtrem(s, t, s);
        if (!fmpz_is_zero(t))
            fmpz_add_ui(s, s, 1UL);

        fmpz_mul(p, p, s);
    }

    fmpz_set(bound, p);
    fmpz_clear(p);
    fmpz_clear(s);
    fmpz_clear(t);
}
