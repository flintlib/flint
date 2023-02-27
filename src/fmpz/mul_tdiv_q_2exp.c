/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

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
fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulong exp)
{
    fmpz c1, c2;
    __mpz_struct * mf;

    c1 = *g;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz_mul_si_tdiv_q_2exp(f, h, c1, exp);
        return;
    }

    c2 = *h;                    /* save h in case it is aliased with f */

    if (c2 == WORD(0))               /* special case, h = 0  */
    {
        fmpz_zero(f);
        return;
    }

    mf = _fmpz_promote(f); /* h is saved, g is already large */

    if (!COEFF_IS_MPZ(c2))      /* g is large, h is small */
        flint_mpz_mul_si(mf, COEFF_TO_PTR(c1), c2);
    else                        /* c1 and c2 are large */
        mpz_mul(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));

    mpz_tdiv_q_2exp(mf, mf, exp);
    _fmpz_demote_val(f);  /* division may make value small */
}
