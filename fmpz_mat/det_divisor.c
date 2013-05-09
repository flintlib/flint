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
#include "fmpq.h"
#include "ulong_extras.h"


void
fmpz_mat_det_divisor(fmpz_t d, const fmpz_mat_t A)
{
    fmpz_mat_t X, B;
    fmpz_t t, u, v, mod;
    len_t i, n;
    int success;

    n = A->r;

    fmpz_mat_init(B, n, 1);
    fmpz_mat_init(X, n, 1);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(mod);

    /* Create a "random" vector */
    for (i = 0; i < n; i++)
    {
        fmpz_set_si(fmpz_mat_entry(B, i, 0), 2*(i % 2) - 1);
    }

    success = fmpz_mat_solve_dixon(X, mod, A, B);

    if (success)
    {
        fmpz_one(d);
        for (i = 0; i < n; i++)
        {
            fmpz_mul(t, d, fmpz_mat_entry(X, i, 0));
            fmpz_fdiv_qr(u, t, t, mod);
            if (!_fmpq_reconstruct_fmpz(u, v, t, mod))
            {
                printf("Exception (fmpz_mat_det_divisor): "
                       "Rational reconstruction failed.\n");
                abort();
            }

            fmpz_mul(d, v, d);
        }
    }
    else
    {
        fmpz_zero(d);
    }

    fmpz_mat_clear(B);
    fmpz_mat_clear(X);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(mod);
}
