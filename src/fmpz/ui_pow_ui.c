/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e)
{
    if (e <= 1)
    {
        fmpz_set_ui(x, e == 0 ? 1 : b);
    }
    else if (e == 2)
    {
        mp_limb_t t[2];
        umul_ppmm(t[1], t[0], b, b);
        fmpz_set_uiui(x, t[1], t[0]);
    }
    else if (b <= 1)
    {
        fmpz_set_ui(x, b);
    }
    else
    {
        ulong bits = FLINT_BIT_COUNT(b);

        if (e * bits <= FLINT_BITS)
        {
            fmpz_set_ui(x, n_pow(b, e));
        }
        else
        {
            __mpz_struct * z = _fmpz_promote(x);
            flint_mpz_set_ui(z, b);
            flint_mpz_pow_ui(z, z, e);
            _fmpz_demote_val(x);
        }
    }
}
