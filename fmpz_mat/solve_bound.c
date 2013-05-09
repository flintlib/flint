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
#include "ulong_extras.h"

void
fmpz_mat_solve_bound(fmpz_t N, fmpz_t D,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    len_t i, j, m, n;
    fmpz_t t, u;

    m = B->r;
    n = B->c;

    fmpz_mat_det_bound(D, A);

    fmpz_init(t);
    fmpz_init(u);

    fmpz_zero(t);

    /* Largest column norm of B */
    for (j = 0; j < n; j++)
    {
        fmpz_zero(u);
        for (i = 0; i < m; i++)
            fmpz_addmul(u, fmpz_mat_entry(B, i, j), fmpz_mat_entry(B, i, j));
        if (fmpz_cmp(t, u) < 0)
            fmpz_set(t, u);
    }

    fmpz_sqrtrem(t, u, t);
    if (!fmpz_is_zero(u))
        fmpz_add_ui(t, t, 1UL);

    fmpz_mul(N, D, t);

    fmpz_clear(t);
    fmpz_clear(u);
}
