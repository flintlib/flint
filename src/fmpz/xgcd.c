/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void
fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
    fmpz_t t1, t2;
    fmpz *f1, *g1;

    if (fmpz_is_zero(f))
    {
        int sign = fmpz_sgn(g);
        fmpz_abs(d, g);
        fmpz_set_ui(a, 0);
        if (sign == 0)
            fmpz_set_ui(b, 0);
        else if (sign > 0)
            fmpz_set_ui(b, 1);
        else
            fmpz_set_si(b, -1);
    }
    else if (fmpz_cmpabs(f, g) == 0)
    {
        if (fmpz_sgn(f) > 0)
        {
            fmpz_set(d, f);
            fmpz_set_ui(a, 1);
        }
        else
        {
            fmpz_neg(d, f);
            fmpz_set_si(a, -1);
        }
        fmpz_set_si(b, 0);
    }
    else
    {
        int sign1 = fmpz_sgn(f);
        int sign2 = fmpz_sgn(g);

        fmpz_init(t1);
        fmpz_init(t2);

        /* support aliasing */
        if (d == f || a == f || sign1 < 0)
        {
            f1 = t1;
            if (sign1 < 0)
                fmpz_neg(f1, f);
            else
                fmpz_set(f1, f);
        }
        else
            f1 = (fmpz *) f;

        if (d == g || a == g || sign2 < 0)
        {
            g1 = t2;
            if (sign2 < 0)
                fmpz_neg(g1, g);
            else
                fmpz_set(g1, g);
        }
        else
            g1 = (fmpz *) g;

        if (fmpz_cmp(f1, g1) < 0)
        {
            fmpz_gcdinv(d, a, f1, g1);
            fmpz_mul(t1, a, f1);
            fmpz_sub(t1, d, t1);
            fmpz_divexact(b, t1, g1);
        }
        else                    /* g < f */
        {
            fmpz_gcdinv(d, b, g1, f1);
            fmpz_mul(t2, b, g1);
            fmpz_sub(t2, d, t2);
            fmpz_divexact(a, t2, f1);
        }

        if (sign1 < 0)
            fmpz_inplace_neg(a);
        if (sign2 < 0)
            fmpz_inplace_neg(b);

        fmpz_clear(t1);
        fmpz_clear(t2);
    }
}

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
