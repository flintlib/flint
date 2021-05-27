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
        slong k, tmp;
        ulong fn = FLINT_ABS(*f);
        ulong gn = FLINT_ABS(*g);

        _fmpz_demote(d);
        _fmpz_demote(a);
        _fmpz_demote(b);

        if (fn == 0 || gn == 0 || fn == gn)
        {
            /* xgcd(0, g) = (|g|, 0, sgn(g)) */
            /* xgcd(f, 0) = (|f|, sgn(f), 0) */
            /* xgcd(±g, g) = (|g|, 0, sgn(g)) */
            *d = (slong) gn + (slong) (fn != gn) * fn;
            *a = (gn == 0) * FLINT_SGN(*f);
            *b = FLINT_SGN(*g);
            return;
        }
        else if (fn == 1 || gn == 1)
        {
            /* xgcd(f, ±1) = (1, 0, ±1)
             * and for g not equal to 0 or ±1, we have
             * xgcd(±1, g) = (1, ±1, 0) */
            *d = 1;
            *a = (fn == 1) * FLINT_SGN(*f);
            *b = (gn == 1) * FLINT_SGN(*g);
            return;
        }
        else if (gn > fn)
        {
            *d = n_xgcd((ulong *) b, (ulong *) a, gn, fn);
        }
        else
            *d = n_xgcd((ulong *) a, (ulong *) b, fn, gn);

        /* We have a solution for +/-(a |f| - b |g|) = d where a, b > 0.
         * Now we want to
         * 1) solve a f + b g = d and
         * 2) make sure |f| =/= 2 d =/= |g|
         * 3) calculate the minimal solution */

        *a *= (2 * (fn >= gn) - 1) * FLINT_SGN(*f);
        *b *= (1 - 2 * (fn >= gn)) * FLINT_SGN(*g); /* we are done with (1) */
        
        if (fn == 2 * *d)
        {
            *b = FLINT_SGN(*g);
            *a = (*d - *b * *g) / *f;
            return;
        }
        else if (gn == 2 * *d)
        {
            *a = FLINT_SGN(*f);
            *b = (*d - *a * *f) / *g;
            return;
        }

        /* The k we want to use lie within d a / g ± 1 / 2 */
        tmp = *g / *d;
        k = *a / tmp;
        k += (FLINT_ABS(*a - (k + 1) * tmp) <= FLINT_ABS(tmp / 2));
        k -= (FLINT_ABS(*a - (k - 1) * tmp) <= FLINT_ABS(tmp / 2));
        *a -= k * tmp;
        *b += k * (*f / *d);

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
