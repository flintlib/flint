/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mpn_extras.h"

int
arf_mul_rnd_down(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec)
{
    mp_size_t xn, yn, zn;
    mp_limb_t hi, lo;
    slong expfix;
    int sgnbit, ret, fix;
    mp_ptr zptr;

    xn = ARF_XSIZE(x);
    yn = ARF_XSIZE(y);
    sgnbit = (xn ^ yn) & 1;
    xn >>= 1;
    yn >>= 1;

    if (yn > xn)
    {
        FLINT_SWAP(arf_srcptr, x, y);
        FLINT_SWAP(mp_size_t, xn, yn);
    }

    /* Either operand is a special value. */
    if (yn == 0)
    {
        arf_mul_special(z, x, y);
        return 0;
    }

    /* xn == yn == 1 */
    if (xn == 1)
    {
        lo = ARF_NOPTR_D(x)[0];
        hi = ARF_NOPTR_D(y)[0];

        umul_ppmm(hi, lo, hi, lo);

        /* Shift so that the top bit is set (branch free version). */
        fix = !(hi >> (FLINT_BITS - 1));
        hi = (hi << fix) | ((lo >> (FLINT_BITS - 1)) & fix);
        lo = (lo << fix);

        ARF_DEMOTE(z);

        if (lo == 0)
        {
            zn = 1;

            if (prec >= FLINT_BITS)
            {
                lo = hi;
                ret = 0;
            }
            else
            {
                lo = MASK_LIMB(hi, FLINT_BITS - prec);
                ret = (lo != hi);
            }
        }
        else
        {
            zn = 2;

            if (prec <= FLINT_BITS)  /* Must be inexact. */
            {
                lo = MASK_LIMB(hi, FLINT_BITS - prec);
                zn = ret = 1;
            }
            else if (prec >= 2 * FLINT_BITS) /* Must be exact. */
            {
                ret = 0;
            }
            else  /* prec < FLINT_BITS < 2 * FLINT_BITS */
            {
                ret = MASK_LIMB(lo, 2 * FLINT_BITS - prec) != lo;
                lo = MASK_LIMB(lo, 2 * FLINT_BITS - prec);

                if (lo == 0)
                {
                    zn = 1;
                    lo = hi;
                }
            }
        }

        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), -fix);
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(zn, sgnbit);

        zptr = ARF_NOPTR_D(z);
        zptr[0] = lo;
        zptr[1] = hi;

        return ret;
    }
    else if (xn == 2)
    {
        mp_limb_t zz[4];
        mp_limb_t x1, x0, y1, y0;

        x0 = ARF_NOPTR_D(x)[0];
        x1 = ARF_NOPTR_D(x)[1];

        if (yn == 2)
        {
            y0 = ARF_NOPTR_D(y)[0];
            y1 = ARF_NOPTR_D(y)[1];

            FLINT_MPN_MUL_2X2(zz[3], zz[2], zz[1], zz[0], x1, x0, y1, y0);

            /* Likely case, must be inexact */
            if (prec <= 2 * FLINT_BITS)
            {
                ARF_DEMOTE(z);

                fix = !(zz[3] >> (FLINT_BITS - 1));
                zz[3] = (zz[3] << fix) | ((zz[2] >> (FLINT_BITS - 1)) & fix);
                zz[2] = (zz[2] << fix) | ((zz[1] >> (FLINT_BITS - 1)) & fix);

                _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), -fix);

                /* Rounding */
                if (prec != 2 * FLINT_BITS)
                {
                    if (prec > FLINT_BITS)
                    {
                        zz[2] &= (LIMB_ONES << (2 * FLINT_BITS - prec));
                    }
                    else if (prec == FLINT_BITS)
                    {
                        zz[2] = 0;
                    }
                    else
                    {
                        zz[3] &= (LIMB_ONES << (FLINT_BITS - prec));
                        zz[2] = 0;
                    }
                }

                if (zz[2] == 0)
                {
                    ARF_XSIZE(z) = ARF_MAKE_XSIZE(1, sgnbit);
                    ARF_NOPTR_D(z)[0] = zz[3];
                }
                else
                {
                    ARF_XSIZE(z) = ARF_MAKE_XSIZE(2, sgnbit);
                    ARF_NOPTR_D(z)[0] = zz[2];
                    ARF_NOPTR_D(z)[1] = zz[3];
                }

                return 1;
            }
        }
        else
        {
            y0 = ARF_NOPTR_D(y)[0];
            FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], x1, x0, y0);
        }

        zn = xn + yn;
        ret = _arf_set_round_mpn(z, &expfix, zz, zn, sgnbit, prec, ARF_RND_DOWN);
        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), expfix);
        return ret;
    }
    else if (yn > MUL_MPFR_MIN_LIMBS && prec != ARF_PREC_EXACT
                && xn + yn > 1.25 * prec / FLINT_BITS
                && xn < MUL_MPFR_MAX_LIMBS)  /* FIXME: proper cutoffs */
    {
        return arf_mul_via_mpfr(z, x, y, prec, ARF_RND_DOWN);
    }
    else
    {
        mp_size_t zn, alloc;
        mp_srcptr xptr, yptr;
        mp_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);
        ARF_GET_MPN_READONLY(yptr, yn, y);

        alloc = zn = xn + yn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);

        ret = _arf_set_round_mpn(z, &expfix, tmp, zn, sgnbit, prec, ARF_RND_DOWN);
        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), expfix);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return ret;
    }
}

int
arf_mul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    slong fix, shift;
    int sgnbit, inexact;

    yn = FLINT_ABS(y->_mp_size);
    xn = ARF_XSIZE(x);
    xn >>= 1;
    sgnbit = ARF_SGNBIT(x) ^ (y->_mp_size < 0);

    /* Either operand is a special value. */
    if (xn == 0 || yn == 0)
    {
        if (arf_is_finite(x))
        {
            arf_zero(z);
        }
        else
        {
            arf_t t;
            arf_init_set_si(t, mpz_sgn(y));
            arf_mul(z, x, t, prec, rnd);
            arf_clear(t);
        }
        return 0;
    }
    else
    {
        mp_size_t zn, alloc;
        mp_srcptr xptr, yptr;
        mp_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);
        yptr = y->_mp_d;

        alloc = zn = xn + yn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);

        shift = yn * FLINT_BITS - (tmp[zn - 1] == 0) * FLINT_BITS;
        zn -= (tmp[zn - 1] == 0);

        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, sgnbit, prec, rnd);

        _fmpz_add_fast(ARF_EXPREF(z), ARF_EXPREF(x), fix + shift);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return inexact;
    }
}

