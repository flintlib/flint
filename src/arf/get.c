/*
    Copyright (C) 2012, 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "double_extras.h"
#include "mpn_extras.h"
#include "fmpq.h"
#include "arf.h"

/* most double: (2^53-1) * 2^971 */
/* least normal: 2^-1022 */
/* least subnormal: 2^-1074 */

static double
huge_double(arf_rnd_t rnd, int negative)
{
    double v;

    if (rnd == ARF_RND_NEAR || arf_rounds_up(rnd, negative))
        v = D_INF;
    else
        v = ldexp(9007199254740991.0, 971);

    return negative ? -v : v;
}

static double
tiny_double(arf_rnd_t rnd, int negative)
{
    double v;

    if (rnd == ARF_RND_NEAR || !arf_rounds_up(rnd, negative))
        v = 0.0;
    else
        v = ldexp(1.0, -1074);

    return negative ? -v : v;
}

double
arf_get_d(const arf_t x, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return 0.0;
        else if (arf_is_pos_inf(x))
            return D_INF;
        else if (arf_is_neg_inf(x))
            return -D_INF;
        else
            return D_NAN;
    }
    else
    {
        arf_t t;
        mp_srcptr tp;
        mp_size_t tn;
        double v;

        /* also catches bignum exponents */
        if (ARF_EXP(x) > 1030 || ARF_EXP(x) < -1080)
        {
            if (fmpz_sgn(ARF_EXPREF(x)) > 0)
                return huge_double(rnd, ARF_SGNBIT(x));
            else
                return tiny_double(rnd, ARF_SGNBIT(x));
        }

        /* allow mpfr to take care of corner cases for now */
        if (ARF_EXP(x) > 1020 || ARF_EXP(x) <= -1020 || rnd == ARF_RND_NEAR)
        {
            mpfr_t xx;
            ARF_GET_MPN_READONLY(tp, tn, x);

            xx->_mpfr_d = (mp_ptr) tp;
            xx->_mpfr_prec = tn * FLINT_BITS;
            xx->_mpfr_sign = ARF_SGNBIT(x) ? -1 : 1;
            xx->_mpfr_exp = ARF_EXP(x);

            return mpfr_get_d(xx, arf_rnd_to_mpfr(rnd));
        }

        arf_init(t);
        arf_set_round(t, x, 53, rnd);
        ARF_GET_MPN_READONLY(tp, tn, t);

        if (tn == 1)
            v = (double)(tp[0]);
        else
            v = (double)(tp[1]) + (double)(tp[0]) * ldexp(1,-32);

        v = d_mul_2exp(v, ARF_EXP(t) - FLINT_BITS);

        if (ARF_SGNBIT(t))
            v = -v;

        arf_clear(t);

        return v;
    }
}

void
arf_get_fmpq(fmpq_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        fmpq_zero(y);
    }
    else if (arf_is_special(x) || !ARF_IS_LAGOM(x))
    {
        flint_throw(FLINT_ERROR, "exception: arf_get_fmpq: cannot convert to rational\n");
    }
    else
    {
        fmpz_t man, exp;
        slong e;

        fmpz_init(man);
        fmpz_init(exp);

        arf_get_fmpz_2exp(man, exp, x);

        e = *exp;

        fmpz_set_ui(fmpq_denref(y), UWORD(1));

        if (e >= 0)
        {
            fmpz_mul_2exp(fmpq_numref(y), man, e);
        }
        else
        {
            fmpz_set(fmpq_numref(y), man);
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), -e);
        }

        fmpz_clear(man);
        fmpz_clear(exp);
    }
}

void
arf_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const arf_t x)
{
    if (arf_is_special(x))
    {
        fmpz_zero(man);
        fmpz_zero(exp);
    }
    else
    {
        mp_srcptr xptr;
        mp_size_t xn;
        int shift;

        ARF_GET_MPN_READONLY(xptr, xn, x);

        shift = flint_ctz(xptr[0]);

        fmpz_sub_ui(exp, ARF_EXPREF(x), xn * FLINT_BITS - shift);

        if (xn == 1)
        {
            if (ARF_SGNBIT(x))
                fmpz_neg_ui(man, xptr[0] >> shift);
            else
                fmpz_set_ui(man, xptr[0] >> shift);
        }
        else
        {
            mpz_ptr z = _fmpz_promote(man);

            if (z->_mp_alloc < xn)
                mpz_realloc(z, xn);

            if (shift == 0)
                flint_mpn_copyi(z->_mp_d, xptr, xn);
            else
                mpn_rshift(z->_mp_d, xptr, xn, shift);

            /* top limb cannot be zero */
            z->_mp_size = ARF_SGNBIT(x) ? -xn : xn;
        }
    }
}

int
arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd)
{
    slong exp;
    int negative, inexact, value, roundup;
    mp_size_t xn, zn;
    mp_srcptr xp;
    mpz_ptr zz;
    mp_ptr zp;
    mp_limb_t v, v2, v3;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            fmpz_zero(z);
            return 0;
        }
        else
        {
            flint_throw(FLINT_ERROR, "arf_get_fmpz: cannot convert infinity or nan to integer\n");
        }
    }

    exp = ARF_EXP(x);
    negative = ARF_SGNBIT(x);

    if (COEFF_IS_MPZ(exp))
    {
        /* tiny */
        if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            if (rnd == ARF_RND_NEAR
                || rnd == ARF_RND_DOWN
                || (rnd == ARF_RND_FLOOR && !negative)
                || (rnd == ARF_RND_CEIL && negative))
            {
                fmpz_zero(z);
            }
            else
            {
                fmpz_set_si(z, negative ? -1 : 1);
            }

            return 1;
        }
        else
        {
            flint_throw(FLINT_ERROR, "arf_get_fmpz: number too large to convert to integer\n");
        }
    }

    /* |x| < 1 */
    if (exp <= 0)
    {
        if (rnd == ARF_RND_NEAR)
        {
            if (exp == 0)
            {
                /* check for the special case +/- 1/2 */
                ARF_GET_MPN_READONLY(xp, xn, x);

                if (xp[xn - 1] < LIMB_TOP || (xn == 1 && xp[xn - 1] == LIMB_TOP))
                    value = 0;
                else
                    value = negative ? -1 : 1;
            }
            else
            {
                value = 0;
            }
        }
        else if (rnd == ARF_RND_DOWN ||
            (rnd == ARF_RND_FLOOR && !negative) ||
            (rnd == ARF_RND_CEIL && negative))
        {
            value = 0;
        }
        else
        {
            value = negative ? -1 : 1;
        }

        _fmpz_demote(z);
        *z = value;
        return 1;
    }

    ARF_GET_MPN_READONLY(xp, xn, x);

    /* Fast case: |x| < 2^31 or 2^63 (must save 1 bit for rounding up!) */
    if (exp < FLINT_BITS)
    {
        v = xp[xn - 1];
        v2 = v >> (FLINT_BITS - exp); /* integral part */
        v3 = v << exp;                /* fractional part (truncated, at least 1 bit) */
        inexact = (xn > 1) || (v3 != 0);

        if (inexact && rnd != ARF_RND_DOWN)
        {
            if (rnd == ARF_RND_NEAR)
            {
                /* round up of fractional part is > 1/2,
                   or if equal to 1/2 and the integral part is odd */
                v2 += (v3 > LIMB_TOP) || (v3 == LIMB_TOP && (xn > 1 || (v2 & 1)));
            }
            else
            {
                v2 += (rnd == ARF_RND_UP) || (negative ^ (rnd == ARF_RND_CEIL));
            }
        }

        if (negative)
            fmpz_neg_ui(z, v2);
        else
            fmpz_set_ui(z, v2);
        return inexact;
    }

    /* |x| >= 1 */

    /* Allocate space for result + 1 extra bit. We need one extra bit
       temporarily to check rounding to nearest. We also need one extra bit
       to round up. */
    zn = (exp + (rnd != ARF_RND_DOWN) + FLINT_BITS - 1) / FLINT_BITS;

    zz = _fmpz_promote(z);
    if (zz->_mp_alloc < zn)
        mpz_realloc2(zz, zn * FLINT_BITS);

    zp = zz->_mp_d;

    if (rnd == ARF_RND_DOWN)
    {
        /* zn is the exact size */
        inexact = _arf_get_integer_mpn(zp, xp, xn, exp);
    }
    else
    {
        zp[zn - 1] = 0;
        inexact = _arf_get_integer_mpn(zp, xp, xn, exp + (rnd == ARF_RND_NEAR));

        if (rnd == ARF_RND_NEAR)
        {
            v = zp[0];
            /* round up if fractional part is >= 1/2 and (there are
               more discarded bits, or the truncated value would be odd) */
            roundup = (v & 1) & (inexact | (v >> 1));
            inexact |= (v & 1);
            mpn_rshift(zp, zp, zn, 1);
            mpn_add_1(zp, zp, zn, roundup);
        }
        else if (inexact && ((rnd == ARF_RND_UP)
            || (negative ^ (rnd == ARF_RND_CEIL))))
        {
            mpn_add_1(zp, zp, zn, 1);
        }

        zn -= (zp[zn - 1] == 0);
    }

    zz->_mp_size = negative ? -zn : zn;
    _fmpz_demote_val(z);
    return inexact;
}

int
arf_get_fmpz_fixed_fmpz(fmpz_t y, const arf_t x, const fmpz_t e)
{
    if (arf_is_special(x))
    {
        return arf_get_fmpz(y, x, ARF_RND_DOWN);
    }
    else
    {
        int inexact;
        fmpz_t exp;
        arf_t tmp;

        fmpz_init(exp);
        fmpz_sub(exp, ARF_EXPREF(x), e);
        arf_init_set_shallow(tmp, x);
        ARF_EXP(tmp) = *exp;
        inexact = arf_get_fmpz(y, tmp, ARF_RND_DOWN);
        fmpz_clear(exp);
        return inexact;
    }
}

int
arf_get_fmpz_fixed_si(fmpz_t y, const arf_t x, slong e)
{
    if (arf_is_special(x))
    {
        return arf_get_fmpz(y, x, ARF_RND_DOWN);
    }
    else
    {
        int inexact;
        fmpz_t exp;
        arf_t tmp;

        fmpz_init(exp);
        fmpz_sub_si(exp, ARF_EXPREF(x), e);
        arf_init_set_shallow(tmp, x);
        ARF_EXP(tmp) = *exp;
        inexact = arf_get_fmpz(y, tmp, ARF_RND_DOWN);
        fmpz_clear(exp);
        return inexact;
    }
}

int
_arf_get_integer_mpn(mp_ptr y, mp_srcptr x, mp_size_t xn, slong exp)
{
    slong bot_exp = exp - xn * FLINT_BITS;

    if (bot_exp >= 0)
    {
        mp_size_t bot_limbs;
        flint_bitcnt_t bot_bits;

        bot_limbs = bot_exp / FLINT_BITS;
        bot_bits = bot_exp % FLINT_BITS;

        flint_mpn_zero(y, bot_limbs);

        if (bot_bits == 0)
            flint_mpn_copyi(y + bot_limbs, x, xn);
        else
            y[bot_limbs + xn] = mpn_lshift(y + bot_limbs, x, xn, bot_bits);

        /* exact */
        return 0;
    }
    else if (exp <= 0)
    {
        /* inexact */
        return 1;
    }
    else
    {
        mp_size_t top_limbs;
        flint_bitcnt_t top_bits;
        mp_limb_t cy;

        top_limbs = exp / FLINT_BITS;
        top_bits = exp % FLINT_BITS;

        if (top_bits == 0)
        {
            flint_mpn_copyi(y, x + xn - top_limbs, top_limbs);
            /* inexact */
            return 1;
        }
        else
        {
            /* can be inexact */
            cy = mpn_rshift(y, x + xn - top_limbs - 1,
                top_limbs + 1, FLINT_BITS - top_bits);

            return (cy != 0) || (top_limbs + 1 != xn);
        }
    }
}

void
arf_get_mag(mag_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        mag_zero(y);
    }
    else if (arf_is_special(x))
    {
        mag_inf(y);
    }
    else
    {
        mp_limb_t t, u;

        ARF_GET_TOP_LIMB(t, x);
        t = (t >> (FLINT_BITS - MAG_BITS)) + LIMB_ONE;

        /* may have rounded up to next power of two */
        u = t >> MAG_BITS;
        /* todo: avoid the addition? check agreement with mag_fast_init_set_arf */
        t = (t >> u) + (u & t);

        _fmpz_add_fast(MAG_EXPREF(y), ARF_EXPREF(x), u);
        MAG_MAN(y) = t;
    }
}

void
arf_get_mag_lower(mag_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        mag_zero(y);
    }
    else if (arf_is_special(x))
    {
        if (arf_is_nan(x))
            mag_zero(y);
        else
            mag_inf(y);
    }
    else
    {
        mp_limb_t t;
        ARF_GET_TOP_LIMB(t, x);
        MAG_MAN(y) = t >> (FLINT_BITS - MAG_BITS);
        _fmpz_set_fast(MAG_EXPREF(y), ARF_EXPREF(x));
    }
}

int
arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd)
{
    int r;

    if (arf_is_special(y))
    {
        if (arf_is_zero(y))
            mpfr_set_zero(x, 0);
        else if (arf_is_pos_inf(y))
            mpfr_set_inf(x, 1);
        else if (arf_is_neg_inf(y))
            mpfr_set_inf(x, -1);
        else
            mpfr_set_nan(x);
        r = 0;
    }
    else if (COEFF_IS_MPZ(*ARF_EXPREF(y)))
    {
        /* Incidentally, COEFF_MIN and COEFF_MAX are exactly the same
           as the minimum and maximum allowed MPFR exponents. We
           assert this to make sure the following code is valid.
           Unfortunately, MPFR provides no convenient function to
           assign a big exponent with automatic underflow/overflow. */
        if (COEFF_MIN > mpfr_get_emin_min() ||
            COEFF_MAX < mpfr_get_emax_max())
        {
            flint_throw(FLINT_ERROR, "unsupported MPFR exponent range: %wd, %wd, %wd, %wd\n",
                COEFF_MIN, mpfr_get_emin_min(), COEFF_MAX, mpfr_get_emax_max());
        }

        if (fmpz_sgn(ARF_EXPREF(y)) > 0)
        {
            if (arf_sgn(y) > 0)
            {
                mpfr_set_inf(x, 1);
                mpfr_nextbelow(x);
            }
            else
            {
                mpfr_set_inf(x, -1);
                mpfr_nextabove(x);
            }

            r = mpfr_mul_2si(x, x, 1, rnd);
        }
        else
        {
            if (arf_sgn(y) > 0)
            {
                mpfr_set_zero(x, 1);
                mpfr_nextabove(x);
            }
            else
            {
                mpfr_set_zero(x, -1);
                mpfr_nextbelow(x);
            }

            r = mpfr_mul_2si(x, x, -1, rnd);
        }
    }
    else
    {
        __mpfr_struct t;
        mp_size_t n;
        mp_srcptr d;

        ARF_GET_MPN_READONLY(d, n, y);

        t._mpfr_d = (mp_ptr) d;
        t._mpfr_exp = ARF_EXP(y);
        t._mpfr_prec = n * FLINT_BITS;
        t._mpfr_sign = ARF_SGNBIT(y) ? -1 : 1;

        r = mpfr_set(x, &t, rnd);
    }

    return r;
}

slong
arf_get_si(const arf_t x, arf_rnd_t rnd)
{
    fmpz_t t;
    slong v;
    fmpz_init(t);
    arf_get_fmpz(t, x, rnd);
    if (!fmpz_fits_si(t))
    {
        flint_throw(FLINT_ERROR, "arf_get_si: result does not fit in a signed slong\n");
    }
    v = fmpz_get_si(t);
    fmpz_clear(t);
    return v;
}
