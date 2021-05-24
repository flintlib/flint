/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

void
fmpz_xgcd_minimal(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
    if (fmpz_is_zero(g))
    {
        slong sgn = fmpz_sgn(f);
        fmpz_abs(d, f);
        fmpz_set_si(a, sgn);
        fmpz_zero(b);
        return;
    }
    else if (fmpz_is_zero(f))
    {
        slong sgn = fmpz_sgn(g);
        fmpz_abs(d, g);
        fmpz_set_si(b, sgn);
        fmpz_zero(a);
        return;
    }

    _fmpz_promote(d);
    _fmpz_promote(a);
    _fmpz_promote(b);

    if (!COEFF_IS_MPZ(*f) && !COEFF_IS_MPZ(*g))  /* both are small */
    {
        ulong tf, tg;
        mpz_t mf, mg;

        tf = fmpz_get_ui(f);
        mf->_mp_d = (mp_limb_t *) &tf;
        mf->_mp_size  = fmpz_sgn(f);

        tg = fmpz_get_ui(g);
        mg->_mp_d = (mp_limb_t *) &tg;
        mg->_mp_size  = fmpz_sgn(g);

        mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), COEFF_TO_PTR(*b),
                   mf, mg);
    }
    else if (!COEFF_IS_MPZ(*f))  /* only f is small */
    {
        ulong tf;
        mpz_t mf;

        tf = fmpz_get_ui(f);
        mf->_mp_d = (mp_limb_t *) &tf;
        mf->_mp_size  = fmpz_sgn(f);

        mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), COEFF_TO_PTR(*b),
                   mf, COEFF_TO_PTR(*g));
    }
    else if (!COEFF_IS_MPZ(*g))  /* only g is small */
    {
        ulong tg;
        mpz_t mg;

        tg = fmpz_get_ui(g);
        mg->_mp_d = (mp_limb_t *) &tg;
        mg->_mp_size  = fmpz_sgn(g);

        mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), COEFF_TO_PTR(*b),
                   COEFF_TO_PTR(*f), mg);
    }
    else /* both are big */
    {
        mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), COEFF_TO_PTR(*b),
                   COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
    }

    _fmpz_demote_val(d);
    _fmpz_demote_val(a);
    _fmpz_demote_val(b);
}
