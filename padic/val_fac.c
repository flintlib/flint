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

#include "padic.h"

ulong padic_val_fac_ui_2(ulong N)
{
    ulong s = 0, t = N;

    do
    {
        t /= 2;
        s += t;
    }
    while (t);

    return s;
}

ulong padic_val_fac_ui(ulong N, const fmpz_t p)
{
    if (fmpz_abs_fits_ui(p))
    {
        ulong q = fmpz_get_ui(p), s = 0, t = N;

        do
        {
            t /= q;
            s += t;
        }
        while (t);

        return s;
    }
    else
    {
        return 0;
    }
}

void padic_val_fac(fmpz_t rop, const fmpz_t op, const fmpz_t p)
{
    fmpz_t t, q, pow;

    if (fmpz_sgn(op) <= 0)
    {
        printf("Exception (padic_val_fac).  op is non-positive.\n");
        abort();
    }

    fmpz_init(t);
    fmpz_init(q);
    fmpz_init(pow);
    fmpz_one(pow);

    do 
    {
        fmpz_mul(pow, pow, p);
        fmpz_fdiv_q(q, op, pow);
        fmpz_add(t, t, q);
    }
    while (!fmpz_is_zero(q));

    fmpz_swap(rop, t);
    fmpz_clear(t);
    fmpz_clear(q);
    fmpz_clear(pow);
}

