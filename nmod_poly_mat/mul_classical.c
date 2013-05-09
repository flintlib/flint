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

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_mul_classical(nmod_poly_mat_t C, const nmod_poly_mat_t A,
    const nmod_poly_mat_t B)
{
    len_t ar, bc, br;
    len_t i, j, k;
    nmod_poly_t t;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (br == 0 || ar == 0 || bc == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, ar, bc, nmod_poly_mat_modulus(A));
        nmod_poly_mat_mul_classical(T, A, B);
        nmod_poly_mat_swap(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    nmod_poly_init(t, nmod_poly_mat_modulus(A));

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            nmod_poly_mul(nmod_poly_mat_entry(C, i, j),
                          nmod_poly_mat_entry(A, i, 0),
                          nmod_poly_mat_entry(B, 0, j));

            for (k = 1; k < br; k++)
            {
                nmod_poly_mul(t, nmod_poly_mat_entry(A, i, k),
                                 nmod_poly_mat_entry(B, k, j));

                nmod_poly_add(nmod_poly_mat_entry(C, i, j),
                              nmod_poly_mat_entry(C, i, j), t);
            }
        }
    }

    nmod_poly_clear(t);
}
