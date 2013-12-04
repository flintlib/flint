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
u_recur(fmpz_t a, fmpz_t b, fmpz_t t, fmpz_t u, fmpz_t v,
        ulong n, const fmpz_t x, int both)
{
    if (n == 0)
    {
        fmpz_one(a);
        if (both)
            fmpz_mul_2exp(b, x, 1);
    }
    else if (n == 1)
    {
        fmpz_mul_2exp(a, x, 1);
        if (both)
        {
            fmpz_mul(b, a, a);
            fmpz_sub_ui(b, b, 1);
        }
    }
    else if (n % 2 == 0)
    {
        u_recur(a, b, t, u, v, (n - 1) / 2, x, 1);

        fmpz_sub(t, b, a);
        fmpz_add(u, b, a);

        /* 2b(xb-a) */
        if (both)
        {
            fmpz_mul(v, x, b);
            fmpz_sub(v, v, a);
            fmpz_mul(v, v, b);
            fmpz_mul_2exp(b, v, 1);
        }

        /* (b+a)(b-a) */
        fmpz_mul(a, t, u);
    }
    else
    {
        u_recur(a, b, t, u, v, (n - 1) / 2, x, 1);

        fmpz_mul(v, x, a);
        fmpz_sub(v, b, v);
        fmpz_mul(v, v, a);

        /* (b+a)*(b-a) */
        if (both)
        {
            fmpz_sub(t, b, a);
            fmpz_add(u, b, a);
            fmpz_mul(b, t, u);
        }

        /* 2a(b-xa) */
        fmpz_mul_2exp(a, v, 1);
    }
}

void
fmpz_chebyshev_u(fmpz_t y, ulong n, const fmpz_t x)
{
    if (fmpz_is_zero(x))
    {
        fmpz_set_si(y, n % 2 ? 0 : (n % 4 ? -1 : 1));
    }
    else if (fmpz_is_pm1(x))
    {
        int sign = (!fmpz_is_one(x) && n % 2);

        fmpz_set_ui(y, n);
        fmpz_add_ui(y, y, 1);
        if (sign)
            fmpz_neg(y, y);
    }
    else
    {
        fmpz_t a, b, t, u, v;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(t);
        fmpz_init(u);
        fmpz_init(v);

        u_recur(a, b, t, u, v, n, x, 0);
        fmpz_swap(y, a);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_clear(v);
    }
}

