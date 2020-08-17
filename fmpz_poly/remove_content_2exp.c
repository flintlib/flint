/*
    Copyright (C) 2016 Elias Tsigaridas
    Copyright (C) 2017 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

slong _fmpz_poly_remove_content_2exp(fmpz * pol, slong len)
{
    slong cont, i, z;

    i = 0;
    while (i < len && fmpz_is_zero(pol + i))
        i++;

    if (i == len)
        return 0;

    cont = fmpz_val2(pol + i);

    for ( ; (i < len) && cont; i++)
    {
        if (!fmpz_is_zero(pol + i))
        {
            z = fmpz_val2(pol + i);
            if (z < cont)
                cont = z;
        }
    }

    if (cont == 0)
        return 0;

    for (i = 0; i < len; i++)
        fmpz_fdiv_q_2exp(pol + i, pol + i, cont);

    return cont;
}
