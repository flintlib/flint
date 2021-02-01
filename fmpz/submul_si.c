/*
    Copyright (C) 2021 Daniel Schultz

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

void fmpz_submul_si(fmpz_t f, const fmpz_t g, slong x)
{
    fmpz F, G;

    G = *g;
    if (x == 0 || G == 0)
        return;

    F = *f;
    if (F == 0)
    {
        fmpz_mul_si(f, g, x);
        fmpz_neg(f, f);
        return;
    }

    if (!COEFF_IS_MPZ(G))
    {
        ulong p1, p0;
        smul_ppmm(p1, p0, G, x);

        if (!COEFF_IS_MPZ(F))
        {
            ulong F1 = FLINT_SIGN_EXT(F);
            sub_ddmmss(p1, p0, F1, F, p1, p0);
            fmpz_set_signed_uiui(f, p1, p0);
        }
        else
        {
            mpz_ptr pF = COEFF_TO_PTR(F);
            sub_ddmmss(p1, p0, UWORD(0), UWORD(0), p1, p0);
            flint_mpz_add_signed_uiui(pF, pF, p1, p0);
        }
    }
    else
    {
        mpz_ptr pG = COEFF_TO_PTR(G);
        mpz_ptr pF = _fmpz_promote_val(f);

        if (x < 0)
            flint_mpz_addmul_ui(pF, pG, -x);
        else
            flint_mpz_submul_ui(pF, pG, x);

        _fmpz_demote_val(f);
    }
}

