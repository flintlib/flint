/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arb.h"
#include "mpn_extras.h"

void
arb_si_pow_ui(arb_t res, slong b, ulong e, slong prec)
{
    arb_ui_pow_ui(res, UI_ABS_SI(b), e, prec);

    if ((e & 1) && b < 0)
        arb_neg(res, res);
}

void
arb_set_round_ui(arb_t res, ulong lo, slong prec)
{
    if (lo == 0)
    {
        arb_zero(res);
    }
    else
    {
        int inexact;

        inexact = _arf_set_round_ui(arb_midref(res), lo, 0, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
        else
            mag_zero(arb_radref(res));
    }
}

void
arb_set_round_uiui(arb_t res, ulong hi, ulong lo, slong prec)
{
    if (hi == 0 && lo == 0)
    {
        arb_zero(res);
    }
    else
    {
        int inexact;
        slong fix;

        inexact = _arf_set_round_uiui(arb_midref(res), &fix, hi, lo, 0, prec, ARB_RND);

        _fmpz_demote(ARF_EXPREF(arb_midref(res)));
        ARF_EXP(arb_midref(res)) = 2 * FLINT_BITS + fix;

        if (inexact)
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
        else
            mag_zero(arb_radref(res));
    }
}

void
arb_ui_pow_ui(arb_t res, ulong a, ulong exp, slong prec)
{
    slong wp, aexp, awidth, trailing, wp_limbs;
    slong exp_fix, alloc, leading;
    int inexact, i, ebits;
    mp_ptr yman, tmp;
    mp_size_t yn;
    mp_limb_t yexp_hi, yexp_lo, aman, aodd, hi, lo;
    ARF_MUL_TMP_DECL

    if (exp <= 2)
    {
        if (exp == 0)
        {
            arb_one(res);
        }
        else if (exp == 1)
        {
            arb_set_round_ui(res, a, prec);
        }
        else
        {
            mp_limb_t hi, lo;
            umul_ppmm(hi, lo, a, a);
            arb_set_round_uiui(res, hi, lo, prec);
        }
        return;
    }

    if (a <= 1)
    {
        arb_set_ui(res, a);
        return;
    }

    aexp = FLINT_BIT_COUNT(a);
    trailing = flint_ctz(a);
    awidth = aexp - trailing;

    /* a = power of two */
    if (awidth == 1)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_ui(t, exp);
        fmpz_mul_ui(t, t, aexp - 1);
        arb_one(res);
        arb_mul_2exp_fmpz(res, res, t);
        fmpz_clear(t);
        return;
    }

    umul_ppmm(hi, lo, awidth, exp);
    if (hi == 0)
        prec = FLINT_MIN(prec, lo);

    wp = prec + 2 * FLINT_BIT_COUNT(exp) + 4;
    wp_limbs = (wp + FLINT_BITS - 1) / FLINT_BITS;

    if (FLINT_BITS == 32 && wp_limbs % 2)
        wp_limbs++;

    /* Algorithm: as long as the result is exact, work with
       powers of the odd part of a as an integer. As soon as we exceed
       wp_limbs, switch to a floating-point product with
       exactly wp_limbs (it is possible that we get trailing
       zero limbs after this point so that the number of limbs
       could be reduced temporarily, but this is not worth
       the trouble). */

    alloc = 4 * (wp_limbs + 2);

    ARF_MUL_TMP_ALLOC(tmp, alloc)
    yman = tmp + 2 * (wp_limbs + 2);

    /* a as a top-aligned mantissa */
    aman = a << (FLINT_BITS - aexp);
    /* a as a bottom-aligned mantissa */
    aodd = a >> trailing;

    /* y = a */
    yn = 1;
    yman[0] = aodd;
    /* yexp will initially be the bottom exponent. we convert
       to a floating-point (top) exponent only when the result
       becomes inexact */
    yexp_lo = trailing;
    /* the exponent can be two words wide (we only need this
       in the inexact case) */
    yexp_hi = 0;

    inexact = 0;
    ebits = FLINT_BIT_COUNT(exp);

    for (i = ebits - 2; i >= 0; i--)
    {
        /* Integer case */
        if (!inexact)
        {
            /* Inline code for small integers */
            if (yn == 1)
            {
                yexp_lo *= 2;
                umul_ppmm(yman[1], yman[0], yman[0], yman[0]);
                yn += (yman[1] != 0);

                if (exp & (UWORD(1) << i))
                {
                    yexp_lo += trailing;

                    if (yn == 1)
                    {
                        umul_ppmm(yman[1], yman[0], yman[0], aodd);
                        yn += (yman[1] != 0);
                    }
                    else
                    {
                        mp_limb_t y0, y1;
                        y0 = yman[0];
                        y1 = yman[1];
                        FLINT_MPN_MUL_2X1(yman[2], yman[1], yman[0], y1, y0, aodd);
                        yn += (yman[2] != 0);
                    }
                }
            }
            else
            {
                /* todo: if 2 * yn is significantly larger than
                   wp_limbs, we might want to go to the floating-point
                   code here */
                yexp_lo *= 2;
                flint_mpn_sqr(tmp, yman, yn);
                yn = 2 * yn;
                yn -= (tmp[yn - 1] == 0);

                if (exp & (UWORD(1) << i))
                {
                    yexp_lo += trailing;
                    yman[yn] = mpn_mul_1(yman, tmp, yn, aodd);
                    yn += (yman[yn] != 0);
                }
                else
                {
                    flint_mpn_copyi(yman, tmp, yn);
                }
            }

            /* convert to floating-point form */
            /* todo: redundant if this is the last iteration */
            if (yn > wp_limbs)
            {
                inexact = 1;
                leading = flint_clz(yman[yn - 1]);
                yexp_lo = yexp_lo + yn * FLINT_BITS - leading;

                if (leading == 0)
                    flint_mpn_copyi(yman, yman + yn - wp_limbs, wp_limbs);
                else
                    mpn_rshift(yman, yman + yn - wp_limbs - 1, wp_limbs + 1, FLINT_BITS - leading);

                yn = wp_limbs;
            }

            continue;
        }

        /* y = y^2: exponent */
        yexp_hi = (yexp_hi << 1) | (yexp_lo >> (FLINT_BITS - 1));
        yexp_lo <<= 1;

        /* special case for 1-limb precision */
        /* note: we must have yn == 1 here if wp_limbs == 1 */
        if (wp_limbs == 1)
        {
            mp_limb_t hi, lo;

            /* y = y^2: mantissa */
            umul_ppmm(hi, lo, yman[0], yman[0]);
            if (!(hi >> (FLINT_BITS - 1)))
            {
                hi = (hi << 1) | (lo >> (FLINT_BITS - 1));
                sub_ddmmss(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, 1);
            }
            yman[0] = hi;

            if (exp & (UWORD(1) << i))
            {
                /* y = y * a: exponent */
                add_ssaaaa(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, aexp);

                /* y = y * a: mantissa */
                umul_ppmm(hi, lo, yman[0], aman);
                if (!(hi >> (FLINT_BITS - 1)))
                {
                    hi = (hi << 1) | (lo >> (FLINT_BITS - 1));
                    sub_ddmmss(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, 1);
                }
                yman[0] = hi;
            }

            continue;
        }

        /* y = y^2: mantissa */
        if (yn >= 25 && yn <= 10000)   /* use mpfr for sqrhi */
        {
            mpfr_t zf, xf;

            zf->_mpfr_d = tmp;
            zf->_mpfr_prec = wp_limbs * FLINT_BITS;
            zf->_mpfr_sign = 1;
            zf->_mpfr_exp = 0;

            xf->_mpfr_d = yman;
            xf->_mpfr_prec = yn * FLINT_BITS;
            xf->_mpfr_sign = 1;
            xf->_mpfr_exp = 0;

            mpfr_sqr(zf, xf, MPFR_RNDD);

            if (zf->_mpfr_exp != 0)
                sub_ddmmss(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, 1);

            yn = wp_limbs;
        }
        else   /* or exactly */
        {
            flint_mpn_sqr(tmp, yman, yn);
            yn = 2 * yn;
        }

        if (yn < wp_limbs)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        /* y = y * a */
        if (exp & (UWORD(1) << i))
        {
            tmp[yn] = mpn_mul_1(tmp + yn - wp_limbs, tmp + yn - wp_limbs, wp_limbs, aman);
            yn++;
            add_ssaaaa(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, aexp);
        }

        /* move result from tmp to yman, and normalize; at this point
           there can be 0, 1 or 2 leading zeros */
        if (tmp[yn - 1] >> (FLINT_BITS - 1))
        {
            flint_mpn_copyi(yman, tmp + yn - wp_limbs, wp_limbs);
        }
        else
        {
            /* assumed so that we can read one extra limb with mpn_rshift */
            /* yn == wp_limbs is only possible if we used mpfr_sqr above
               and there was no multiplication by a, but in that case
               there are 0 leading zeros anyway */
            if (yn == wp_limbs)
                flint_throw(FLINT_ERROR, "(%s)\n", __func__);

            if (tmp[yn - 1] >> (FLINT_BITS - 2))
            {
                mpn_rshift(yman, tmp + yn - wp_limbs - 1, wp_limbs + 1, FLINT_BITS - 1);
                sub_ddmmss(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, 1);
            }
            else
            {
                mpn_rshift(yman, tmp + yn - wp_limbs - 1, wp_limbs + 1, FLINT_BITS - 2);
                sub_ddmmss(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, 2);
            }
        }

        yn = wp_limbs;
    }

    /* convert bottom exponent to floating-point (top) exponent */
    if (!inexact)
        yexp_lo = yexp_lo + yn * FLINT_BITS;

    inexact |= _arf_set_round_mpn(arb_midref(res), &exp_fix, yman, yn, 0, prec, ARF_RND_DOWN);
    if (exp_fix)
        sub_ddmmss(yexp_hi, yexp_lo, yexp_hi, yexp_lo, 0, -exp_fix);
    fmpz_set_uiui(ARF_EXPREF(arb_midref(res)), yexp_hi, yexp_lo);

    if (inexact)
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec - 1);
    else
        mag_zero(arb_radref(res));

    ARF_MUL_TMP_FREE(tmp, alloc)
}
