/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

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

void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        if (d >= 0)
            d = d >> FLINT_MIN(exp, FLINT_BITS - 2);
        else
            d = -((-d) >> FLINT_MIN(exp, FLINT_BITS - 2));

        fmpz_set_si(f, d);
    }
    else  /*g is large */
    {
        __mpz_struct * mpz_ptr = _fmpz_promote(f);  /* g is already large */
        mpz_tdiv_q_2exp(mpz_ptr, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}
