/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_scalar_submul_si_2exp(fmpz * vec1, const fmpz * vec2, slong len2,
                                slong c, ulong exp)
{
    slong i;
    fmpz_t temp;

    if (c == 0)
        return;                 /* nothing to add */

    if (exp == 0)               /* just do submul */
    {
        _fmpz_vec_scalar_submul_si(vec1, vec2, len2, c);
        return;
    }

    fmpz_init(temp);

    if (c == 1)                 /* scalar is 1, just subtract c * 2^exp times c */
    {
        for (i = 0; i < len2; i++)
        {
            fmpz_mul_2exp(temp, vec2 + i, exp);
            fmpz_sub(vec1 + i, vec1 + i, temp);
        }
    }
    else if (c == -1)           /* scalar is -1, add c * 2^exp */
    {
        for (i = 0; i < len2; i++)
        {
            fmpz_mul_2exp(temp, vec2 + i, exp);
            fmpz_add(vec1 + i, vec1 + i, temp);
        }
    }
    else                        /* generic case */
    {
        if (c >= 0)
        {
            for (i = 0; i < len2; i++)
            {
                fmpz_mul_2exp(temp, vec2 + i, exp);
                fmpz_submul_ui(vec1 + i, temp, c);
            }
        }
        else
        {
            for (i = 0; i < len2; i++)
            {
                fmpz_mul_2exp(temp, vec2 + i, exp);
                fmpz_addmul_ui(vec1 + i, temp, -c);
            }
        }
    }

    fmpz_clear(temp);
}
