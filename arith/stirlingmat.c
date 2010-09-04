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

#include "flint.h"
#include "fmpz.h"
#include "arith.h"


void
_fmpz_stirling_mat(fmpz ** rows, long nn, int kind)
{
    long n, k;

    fmpz * row;
    fmpz * prev_row;

    for (n = 0; n < nn; n++)
    {
        prev_row = row;
        row = rows[n];

        /* Leftmost column and diagonal */
        fmpz_set_ui(row, n == 0);
        fmpz_set_ui(row + n, 1UL);

        for (k = 1; k < n; k++)
        {
            fmpz_set(row+k, prev_row+(k-1));
            switch (kind)
            {
            case 0:
                fmpz_addmul_ui(row+k, prev_row+k, n - 1UL);
                break;
            case 1:
                fmpz_submul_ui(row+k, prev_row+k, n - 1UL);
                break;
            case 2:
                fmpz_addmul_ui(row+k, prev_row+k, k);
                break;
            }
        }
    }
}

void
fmpz_stirling1u_mat(fmpz ** rows, long n)
{
    _fmpz_stirling_mat(rows, n, 0);
}

void
fmpz_stirling1_mat(fmpz ** rows, long n)
{
    _fmpz_stirling_mat(rows, n, 1);
}

void
fmpz_stirling2_mat(fmpz ** rows, long n)
{
    _fmpz_stirling_mat(rows, n, 2);
}
