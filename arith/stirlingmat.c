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
_fmpz_stirling_next_row(fmpz * new, fmpz * prev, long n, long klen, int kind)
{
    long k;
    fmpz_t t, u;

    if (n == 0)
    {
        fmpz_set_ui(new, 1UL);
        return;
    }

    if (klen < 1)
        return;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_set_ui(new, 0UL);

    if (klen > n)
        fmpz_set_ui(new+n, 1UL);

    for (k = 1; k < FLINT_MIN(n, klen); k++)
    {
        fmpz_set(u, prev+k);
        fmpz_set(new + k, t);
        switch (kind)
        {
        case 0:
            fmpz_addmul_ui(new+k, u, n-1UL);
            break;
        case 1:
            fmpz_submul_ui(new+k, u, n - 1UL);
            break;
        case 2:
            fmpz_addmul_ui(new+k, u, k);
            break;
        }
        fmpz_set(t, u);
    }
    fmpz_clear(t);
    fmpz_clear(u);
}

void
_fmpz_stirling_mat(fmpz ** rows, long nn, int kind)
{
    long n;

    if (nn > 0)
        fmpz_set_ui(rows[0], 1UL);
    for (n = 1; n < nn; n++)
        _fmpz_stirling_next_row(rows[n], rows[n-1], n, n+1, kind);
}

void fmpz_stirling1u_vec_next(fmpz * row, fmpz * prev, long n, long klen)
{
    _fmpz_stirling_next_row(row, prev, n, klen, 0);
}

void fmpz_stirling1_vec_next(fmpz * row, fmpz * prev, long n, long klen)
{
    _fmpz_stirling_next_row(row, prev, n, klen, 1);
}

void fmpz_stirling2_vec_next(fmpz * row, fmpz * prev, long n, long klen)
{
    _fmpz_stirling_next_row(row, prev, n, klen, 2);
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
