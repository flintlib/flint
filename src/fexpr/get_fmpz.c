/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_extras.h"
#include "fexpr.h"

int
fexpr_get_fmpz(fmpz_t c, const fexpr_t x)
{
    ulong head = x->data[0];

    if (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_INT)
    {
        _fmpz_demote(c);
        *c = ((slong) head) >> FEXPR_TYPE_BITS;
    }
    else
    {
        slong nlimbs;
        int negative;

        nlimbs = FEXPR_SIZE(head) - 1;

        if (FEXPR_TYPE(head) == FEXPR_TYPE_BIG_INT_POS)
        {
            negative = 0;
        }
        else if (FEXPR_TYPE(head) == FEXPR_TYPE_BIG_INT_NEG)
        {
            negative = 1;
        }
        else
        {
            return 0;
        }

        if (nlimbs == 1 && x->data[1] <= COEFF_MAX)
        {
            _fmpz_demote(c);
            *c = negative ? (-(slong) x->data[1]) : x->data[1];
        }
        else
        {
            fmpz_set_mpn_large(c, x->data + 1, nlimbs, negative);
        }
    }

    return 1;
}
