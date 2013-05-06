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
fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X,
                                    const fmpz_mat_t Xmod, const fmpz_t mod)
{
    fmpz_t num, den, t, u, d;
    long i, j;

    int success = 1;

    fmpz_init(num);
    fmpz_init(den);
    fmpz_init(d);
    fmpz_init(t);
    fmpz_init(u);

    fmpz_one(d);

    for (i = 0; i < Xmod->r; i++)
    {
        for (j = 0; j < Xmod->c; j++)
        {
            /* TODO: handle various special cases efficiently; zeros,
                     small integers, etc. */
            fmpz_mul(t, d, fmpz_mat_entry(Xmod, i, j));
            fmpz_fdiv_qr(u, t, t, mod);

            success = _fmpq_reconstruct_fmpz(num, den, t, mod);

            fmpz_mul(den, den, d);
            fmpz_set(d, den);

            if (!success)
                goto cleanup;

            fmpz_set(fmpq_mat_entry_num(X, i, j), num);
            fmpz_set(fmpq_mat_entry_den(X, i, j), den);
            fmpq_canonicalise(fmpq_mat_entry(X, i, j));
        }
    }

cleanup:
    fmpz_clear(num);
    fmpz_clear(den);
    fmpz_clear(d);
    fmpz_clear(t);
    fmpz_clear(u);

    return success;
}
