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
#include "fmpq.h"
#include "fmpq_mat.h"


void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
{
    long i, j, k;

    if (A == C || B == C)
    {
        printf("Exception (fmpq_mat_mul_direct). Aliasing not implemented.\n");
        abort();
    }

    if (A->c == 0)
    {
        fmpq_mat_zero(C);
        return;
    }

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            fmpq_mul(fmpq_mat_entry(C, i, j),
                     fmpq_mat_entry(A, i, 0),
                     fmpq_mat_entry(B, 0, j));

            for (k = 1; k < A->c; k++)
            {
                fmpq_addmul(fmpq_mat_entry(C, i, j),
                            fmpq_mat_entry(A, i, k),
                            fmpq_mat_entry(B, k, j));
            }
        }
    }
}
