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

    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "flint.h"
#include "mpfr_mat.h"

void
mpfr_mat_mul_classical(mpfr_mat_t C, const mpfr_mat_t A, const mpfr_mat_t B,
                       mpfr_rnd_t rnd)
{
    slong ar, bc, br;
    slong i, j, k;
    mpfr_t tmp;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (C == A || C == B)
    {
        mpfr_mat_t t;
        mpfr_mat_init(t, ar, bc, C->prec);
        mpfr_mat_mul_classical(t, A, B, rnd);
        mpfr_mat_swap(C, t);
        mpfr_mat_clear(t);
        return;
    }

    if (C->r != ar || C->c != bc)
    {
        flint_printf
            ("Exception (mpfr_mat_mul_classical). Incompatible dimensions.\n");
        abort();
    }

    if (br == 0)
    {
        mpfr_mat_zero(C);
        return;
    }

    mpfr_init2(tmp, C->prec);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mpfr_mul(mpfr_mat_entry(C, i, j), mpfr_mat_entry(A, i, 0),
                     mpfr_mat_entry(B, 0, j), rnd);

            for (k = 1; k < br; k++)
            {
                mpfr_mul(tmp, mpfr_mat_entry(A, i, k), mpfr_mat_entry(B, k, j),
                         rnd);
                mpfr_add(mpfr_mat_entry(C, i, j), mpfr_mat_entry(C, i, j), tmp,
                         rnd);
            }
        }
    }

    mpfr_clear(tmp);
}
