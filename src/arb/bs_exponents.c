/*
    Copyright (C) 2014, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb/impl.h"

/* The exponent table of the binary splittings with a power table (the
   exponential series of mp_real, padic_radix, arb_hypgeom/gamma_fmpq). */

slong _arb_compute_bs_exponents(slong * tab, slong n)
{
    slong a, b, aa, /* ab, */ ba, bb, length;

    if (n == 1)
    {
        tab[0] = 1;
        return 1;
    }

    if (n == 2 || n == 3 || n == 4)
    {
        tab[0] = 1;
        tab[1] = 2;
        return 2;
    }

    if (n == 6)
    {
        tab[0] = 1;
        tab[1] = 2;
        tab[2] = 3;
        return 3;
    }

    /* first binary splitting call */
    a = n >> 1;
    b = n - (n >> 1);
    tab[0] = a;
    length = 1;

    for (;;)
    {
        /* split a -> aa, ab  and b -> ba, bb */
        aa = a >> 1;
        /* ab = a - aa; */
        ba = b >> 1;
        bb = b - ba;

        tab[length] = ba;
        length++;

        /* at length 3, we split into 2, 1 (and maybe also 2, 2) */
        if (ba == 3)
        {
            tab[length] = 2;
            tab[length + 1] = 1;
            length += 2;
            break;
        }

        /* stop if we have reached 1, or if at 2 and the
           length is a power of 2 (in which case we never reach 1) */
        if (ba == 1 || (ba == 2 && (n & (n-1)) == 0))
            break;

        /* if left and right lengths are different, also add the
           right length */
        if (aa != ba && aa != 1)
        {
            tab[length] = aa;
            length++;
        }

        a = aa;
        b = bb;
    }

    /* we always include x^1 in the table, even if the binary splitting
       terminates at step length 2 */
    if (tab[length-1] != 1)
    {
        tab[length] = 1;
        length++;
    }

    /* reverse table */
    for (a = 0; a < length / 2; a++)
    {
        b = tab[a];
        tab[a] = tab[length - a - 1];
        tab[length - a - 1] = b;
    }

    return length;
}

/* just do a linear search */
slong _arb_get_exp_pos(const slong * tab, slong step)
{
    slong i;

    for (i = 0; ; i++)
    {
        if (tab[i] == step)
            return i;

        if (tab[i] == 0)
        {
            flint_throw(FLINT_ERROR, "ERROR: exponent %wd not in table!\n", step);
        }
    }
}
