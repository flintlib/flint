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
#include "nmod_mat.h"


void
fmpz_mat_CRT_ui(fmpz_mat_t res, const fmpz_mat_t mat1,
                        const fmpz_t m1, const nmod_mat_t mat2, int sign)
{
    len_t i, j;
    mp_limb_t c;
    mp_limb_t m2 = mat2->mod.n;
    mp_limb_t m2inv = mat2->mod.ninv;
    fmpz_t m1m2;

    c = fmpz_fdiv_ui(m1, m2);
    c = n_invmod(c, m2);

    if (c == 0)
    {
        printf("Exception (fmpz_mat_CRT_ui). m1 not invertible modulo m2.\n");
        abort();
    }

    fmpz_init(m1m2);
    fmpz_mul_ui(m1m2, m1, m2);

    for (i = 0; i < mat1->r; i++)
    {
        for (j = 0; j < mat1->c; j++)
            _fmpz_CRT_ui_precomp(fmpz_mat_entry(res, i, j),
                                 fmpz_mat_entry(mat1, i, j), m1,
                    nmod_mat_entry(mat2, i, j), m2, m2inv, m1m2, c, sign);
    }

    fmpz_clear(m1m2);
}
