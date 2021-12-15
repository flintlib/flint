/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 < WORD(0))
            fmpz_addmul_ui(f, h, -c1);
        else
            fmpz_submul_ui(f, h, c1);
    }
    else
    {
        fmpz c2 = *h;
        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            if (c2 < WORD(0))
                fmpz_addmul_ui(f, g, -c2);
            else
                fmpz_submul_ui(f, g, c2);
        }
        else                    /* both g and h are large */
        {
            __mpz_struct * mf = _fmpz_promote_val(f);
            mpz_submul(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);    /* cancellation may have occurred */
        }
    }
}
