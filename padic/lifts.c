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

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include "fmpz.h"

len_t * _padic_lifts_exps(len_t *n, len_t N)
{
    len_t *a, i;

    *n = FLINT_CLOG2(N) + 1;

    a = flint_malloc((*n) * sizeof(len_t));
    for (a[i = 0] = N; a[i] > 1; i++)
        a[i + 1] = (a[i] + 1) / 2;

    return a;
}

void _padic_lifts_pows(fmpz *pow, const len_t *a, len_t n, const fmpz_t p)
{
    if (n == 1)
    {
        fmpz_set(pow + 0, p);
    }
    else  /* n > 1 */
    {
        len_t i = n - 1;
        fmpz_t t = {1L};

        fmpz_set(pow + i, p);

        for (i--; i >= 1; i--)
        {
            if (a[i] & 1L)
            {
                fmpz_mul(pow + i, t, pow + (i + 1));
                fmpz_mul(t, t, t);
            }
            else
            {
                fmpz_mul(t, t, pow + (i + 1));
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
            }
        }
        {
            if (a[i] & 1L)
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }
        fmpz_clear(t);
    }
}

