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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"


void _fmpz_mat_solve_fflu_precomp(fmpz * b, fmpz ** const LU, long n)
{
    long i, j;

    /* Fraction-free forward substitution */
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            fmpz_mul(&b[j], &b[j], &LU[i][i]);
            fmpz_submul(&b[j], &LU[j][i], &b[i]);
            if (i > 0)
                fmpz_divexact(&b[j], &b[j], &LU[i-1][i-1]);
        }
    }

    /* Fraction-free backward substitution */
    for (i = n - 2; i >= 0; i--)
    {
        fmpz_mul(&b[i], &b[i], &LU[n-1][n-1]);
        for (j = i + 1; j < n; j++)
            fmpz_submul(&b[i], &b[j], &LU[i][j]);
        fmpz_divexact(&b[i], &b[i], &LU[i][i]);
    }
}
