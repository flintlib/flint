/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void
fmpz_xgcd_canonical_bezout(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
    /* check aliasing */
    if (d == f || a == f || b == f || d == g || a == g || b == g)
    {
        fmpz_t d2, a2, b2;
        fmpz_init(d2);
        fmpz_init(a2);
        fmpz_init(b2);
        fmpz_xgcd_canonical_bezout(d2, a2, b2, f, g);
        fmpz_swap(d, d2);
        fmpz_swap(a, a2);
        fmpz_swap(b, b2);
        fmpz_clear(d2);
        fmpz_clear(a2);
        fmpz_clear(b2);
        return;
    }


    if (!COEFF_IS_MPZ(*f) && !COEFF_IS_MPZ(*g))  /* both are small */
    {
        ulong fn = FLINT_ABS(*f);
        ulong gn = FLINT_ABS(*g);

        _fmpz_demote(d);
        _fmpz_demote(a);
        _fmpz_demote(b);

        if (fn == 0 || gn == 0)
        {
            /* xgcd(0, g) = (|g|, 0, sgn(g)) */
            /* xgcd(f, 0) = (|f|, sgn(f), 0) */
            *d = (slong) gn + (slong) (fn != gn) * fn;
            *a = (gn == 0) * FLINT_SGN(*f);
            *b = FLINT_SGN(*g);
            return;
        }

        *d = mpn_gcdext_1(a, b, fn, gn);

        *a *= FLINT_SGN(*f);
        *b *= FLINT_SGN(*g);
        return;
    }
    else if (!COEFF_IS_MPZ(*f))  /* only f is small */
    {
        mpz_t mf;
        ulong tf = FLINT_ABS(*f);

        mf->_mp_d = (mp_limb_t *) &tf;
        mf->_mp_size  = fmpz_sgn(f);

        _fmpz_promote(d);
        _fmpz_promote(a);
        _fmpz_promote(b);

        mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), COEFF_TO_PTR(*b),
                   mf, COEFF_TO_PTR(*g));
    }
    else if (!COEFF_IS_MPZ(*g))  /* only g is small */
    {
        mpz_t mg;
        ulong tg = FLINT_ABS(*g);

        mg->_mp_d = (mp_limb_t *) &tg;
        mg->_mp_size  = fmpz_sgn(g);

        _fmpz_promote(d);
        _fmpz_promote(a);
        _fmpz_promote(b);

        mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), COEFF_TO_PTR(*b),
                   COEFF_TO_PTR(*f), mg);
    }
    else /* both are big */
    {
        _fmpz_promote(d);
        _fmpz_promote(a);
        _fmpz_promote(b);

        mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), COEFF_TO_PTR(*b),
                   COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
    }

    _fmpz_demote_val(d);
    _fmpz_demote_val(a);
    _fmpz_demote_val(b);
}
