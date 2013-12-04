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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpz.h"

static void
t_recur(fmpz_t a, fmpz_t b, ulong n, const fmpz_t x, int both)
{
    if (n == 0)
    {
        fmpz_one(a);
        if (both)
            fmpz_set(b, x);
    }
    else if (n == 1)
    {
        fmpz_set(a, x);
        if (both)
            fmpz_one(b);
    }
    else if (n % 2 == 0)
    {
        t_recur(a, b, n / 2, x, both);

        if (both)
        {
            fmpz_mul(b, b, a);
            fmpz_mul_2exp(b, b, 1);
            fmpz_sub(b, b, x);
        }

        fmpz_mul(a, a, a);
        fmpz_mul_2exp(a, a, 1);
        fmpz_sub_ui(a, a, 1);
    }
    else
    {
        t_recur(a, b, n / 2 + 1, x, 1);

        fmpz_mul(a, a, b);
        fmpz_mul_2exp(a, a, 1);
        fmpz_sub(a, a, x);

        if (both)
        {
            fmpz_mul(b, b, b);
            fmpz_mul_2exp(b, b, 1);
            fmpz_sub_ui(b, b, 1);
        }
    }
}

void
fmpz_chebyshev_t(fmpz_t y, ulong n, const fmpz_t x)
{
    if (fmpz_is_zero(x))
    {
        fmpz_set_si(y, n % 2 ? 0 : (n % 4 ? -1 : 1));
    }
    else if (fmpz_is_pm1(x))
    {
        if (fmpz_is_one(x))
            fmpz_one(y);
        else
            fmpz_set_si(y, n % 2 ? -1 : 1);
    }
    else
    {
        fmpz_t a, b;

        fmpz_init(a);
        fmpz_init(b);

        t_recur(a, b, n, x, 0);
        fmpz_swap(y, a);

        fmpz_clear(a);
        fmpz_clear(b);
    }
}

