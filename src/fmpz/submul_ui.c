/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz F, G;

    G = *g;
    if (x == 0 || G == 0)
        return;

    F = *f;
    if (F == 0)
    {
        fmpz_mul_ui(f, g, x);
        fmpz_neg(f, f);
        return;
    }

    if (!COEFF_IS_MPZ(G))
    {
        ulong p1, p0;

        if (x <= WORD_MAX)
        {
            smul_ppmm(p1, p0, -G, x);
        }
        else
        {
            umul_ppmm(p1, p0, FLINT_ABS(G), x);
            if (G > 0)
            {
                p1 = -p1 - (p0 != 0);
                p0 = -p0;
            }
        }

        if (!COEFF_IS_MPZ(F))
        {
            ulong F1 = FLINT_SIGN_EXT(F);
            add_ssaaaa(p1, p0, p1, p0, F1, F);
            fmpz_set_signed_uiui(f, p1, p0);
        }
        else
        {
            mpz_ptr pF = COEFF_TO_PTR(F);
            flint_mpz_add_signed_uiui(pF, pF, p1, p0);
            _fmpz_demote_val(f);
        }
    }
    else
    {
        mpz_ptr pG = COEFF_TO_PTR(G);
        mpz_ptr pF = _fmpz_promote_val(f);

        flint_mpz_submul_ui(pF, pG, x);
        _fmpz_demote_val(f);
    }
}
