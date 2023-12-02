/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

ulong padic_val_fac_ui(ulong N, const fmpz_t prime)
{
    if (fmpz_abs_fits_ui(prime))
    {
        const ulong p = fmpz_get_ui(prime);

        ulong s = 0, t = N;

        do
        {
            t /= p;
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
    fmpz_t s, t;

    if (fmpz_sgn(op) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (padic_val_fac).  op is negative.\n");
    }

    fmpz_init(s);
    fmpz_init_set(t, op);

    do
    {
        fmpz_fdiv_q(t, t, p);
        fmpz_add(s, s, t);
    }
    while (!fmpz_is_zero(t));

    fmpz_swap(rop, s);

    fmpz_clear(s);
    fmpz_clear(t);
}

