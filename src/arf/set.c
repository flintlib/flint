/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "double_extras.h"
#include "arf.h"

void
arf_set(arf_t y, const arf_t x)
{
    if (x != y)
    {
        /* Fast path */
        if (!COEFF_IS_MPZ(ARF_EXP(x)) && !COEFF_IS_MPZ(ARF_EXP(y)))
            ARF_EXP(y) = ARF_EXP(x);
        else
            fmpz_set(ARF_EXPREF(y), ARF_EXPREF(x));

        /* Fast path */
        if (!ARF_HAS_PTR(x))
        {
            ARF_DEMOTE(y);
            (y)->d = (x)->d;
        }
        else
        {
            mp_ptr yptr;
            mp_srcptr xptr;
            mp_size_t n;

            ARF_GET_MPN_READONLY(xptr, n, x);
            ARF_GET_MPN_WRITE(yptr, n, y);
            flint_mpn_copyi(yptr, xptr, n);
        }

        /* Important. */
        ARF_XSIZE(y) = ARF_XSIZE(x);
    }
}

void
arf_set_d(arf_t x, double v)
{
#if FLINT_BITS == 64
    mp_limb_t h, sign, exp, frac;
    slong real_exp, real_man;
    union { double uf; mp_limb_t ul; } u;

    u.uf = v;
    h = u.ul;
    sign = h >> 63;
    exp = (h << 1) >> 53;
    frac = (h << 12) >> 12;

    if (exp == 0 && frac == 0)
    {
        arf_zero(x);
        return;
    }
    else if (exp == 0x7ff)
    {
        if (frac == 0)
        {
            if (sign)
                arf_neg_inf(x);
            else
                arf_pos_inf(x);
        }
        else
        {
            arf_nan(x);
        }
        return;
    }

    /* handle subnormals */
    if (exp == 0 && frac != 0)
    {
        int exp2;
        v = frexp(v, &exp2);
        u.uf = v;
        h = u.ul;
        sign = h >> 63;
        exp = (h << 1) >> 53;
        frac = (h << 12) >> 12;
        exp += exp2;
    }

    real_exp = exp - 1023 - 52;

    frac |= (UWORD(1) << 52);
    real_man = sign ? (-frac) : frac;

    arf_set_si_2exp_si(x, real_man, real_exp);
#else
    mpfr_t t;
    mp_limb_t tmp[2];

    t->_mpfr_prec = 53;
    t->_mpfr_sign = 1;
    t->_mpfr_exp = 0;
    t->_mpfr_d = tmp;

    mpfr_set_d(t, v, MPFR_RNDD);

    arf_set_mpfr(x, t);
#endif
}

void
arf_set_mpfr(arf_t x, const mpfr_t y)
{
    if (!mpfr_regular_p(y))
    {
        if (mpfr_zero_p(y))
            arf_zero(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) > 0)
            arf_pos_inf(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) < 0)
            arf_neg_inf(x);
        else
            arf_nan(x);
    }
    else
    {
        mp_size_t n = (y->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS;
        arf_set_mpn(x, y->_mpfr_d, n, y->_mpfr_sign < 0);
        fmpz_set_si(ARF_EXPREF(x), y->_mpfr_exp);
    }
}

void
arf_set_mpn(arf_t y, mp_srcptr x, mp_size_t xn, int sgnbit)
{
    unsigned int leading;
    mp_size_t yn, xn1;
    mp_ptr yptr;
    mp_limb_t bot;

    xn1 = xn;

    while (x[0] == 0)
    {
        x++;
        xn--;
    }

    leading = flint_clz(x[xn - 1]);

    bot = x[0];

    /* This works when leading == 0, since x[0] != 0. */
    yn = xn - ((bot << leading) == 0);

    ARF_GET_MPN_WRITE(yptr, yn, y);
    ARF_XSIZE(y) |= sgnbit;

    if (leading == 0)
    {
        flint_mpn_copyi(yptr, x, xn);
    }
    else if (xn == yn)
    {
        mpn_lshift(yptr, x, yn, leading);
    }
    else
    {
        mpn_lshift(yptr, x + 1, yn, leading);
        yptr[0] |= (bot >> (FLINT_BITS - leading));
    }

    fmpz_set_ui(ARF_EXPREF(y), xn1 * FLINT_BITS - leading);
}

int
_arf_set_mpn_fixed(arf_t z, mp_srcptr xp, mp_size_t xn,
        mp_size_t fixn, int negative, slong prec, arf_rnd_t rnd)
{
    slong exp, exp_shift;
    int inexact;

    exp = (slong)(xn - fixn) * FLINT_BITS;

    while (xn > 0 && xp[xn-1] == 0)
    {
        xn--;
        exp -= FLINT_BITS;
    }

    if (xn == 0)
    {
        arf_zero(z);
        return 0;
    }
    else
    {
        inexact = _arf_set_round_mpn(z, &exp_shift, xp, xn, negative, prec, rnd);
        fmpz_set_si(ARF_EXPREF(z), exp + exp_shift);
        return inexact;
    }
}
