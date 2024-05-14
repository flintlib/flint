/*
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Arb authors
    Copyright (C) 2019 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_vec.h"

void
_arb_vec_printn(arb_srcptr vec, slong len, slong ndigits, ulong flags)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        arb_printn(vec + i, ndigits, flags);
        if (i < len - 1)
            flint_printf(", ");
    }
}

void
_arb_vec_printd(arb_srcptr vec, slong len, slong ndigits)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        arb_printd(vec + i, ndigits);
        if (i < len - 1)
        {
            flint_printf(", ");
        }
    }
    flint_printf("\n");
}
