/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

FLINT_TLS_PREFIX mp_ptr __arf_add_tmp = NULL;
FLINT_TLS_PREFIX slong __arf_add_alloc = 0;

void _arf_add_tmp_cleanup(void)
{
    flint_free(__arf_add_tmp);
    __arf_add_tmp = NULL;
    __arf_add_alloc = 0;
}

/* Assumptions: top limbs of x and y nonzero. */
int
_arf_add_mpn(arf_t z, mp_srcptr xp, mp_size_t xn, int xsgnbit, const fmpz_t xexp,
                      mp_srcptr yp, mp_size_t yn, int ysgnbit, flint_bitcnt_t shift,
                      slong prec, arf_rnd_t rnd)
{
    mp_size_t wn, zn, zn_original, alloc, xbase, wbase;
    mp_size_t shift_limbs;
    flint_bitcnt_t shift_bits;
    int inexact;
    slong fix;
    mp_limb_t cy;
    mp_ptr tmp;
    ARF_ADD_TMP_DECL

    /* very fast case */
    if (xn == 1 && yn == 1 && shift < FLINT_BITS - 1)
    {
        mp_limb_t hi, lo, xhi, xlo, yhi, ylo;

        xhi = xp[0];
        yhi = yp[0];

        xlo = xhi << (FLINT_BITS - 1);
        xhi = xhi >> 1;

        ylo = yhi << (FLINT_BITS - (shift + 1));
        yhi = yhi >> (shift + 1);

        if (xsgnbit == ysgnbit)
        {
            add_ssaaaa(hi, lo, xhi, xlo, yhi, ylo);
        }
        else
        {
            if (xhi > yhi)
            {
                sub_ddmmss(hi, lo, xhi, xlo, yhi, ylo);
            }
            else if (xhi < yhi)
            {
                sub_ddmmss(hi, lo, yhi, ylo, xhi, xlo);
                xsgnbit = ysgnbit;
            }
            else
            {
                if (xlo > ylo)
                {
                    sub_ddmmss(hi, lo, xhi, xlo, yhi, ylo);
                }
                else if (xlo < ylo)
                {
                    sub_ddmmss(hi, lo, yhi, ylo, xhi, xlo);
                    xsgnbit = ysgnbit;
                }
                else
                {
                    arf_zero(z);
                    return 0;
                }
            }
        }

        inexact = _arf_set_round_uiui(z, &fix, hi, lo, xsgnbit, prec, rnd);
        _fmpz_add_fast(ARF_EXPREF(z), xexp, fix + 1);
        return inexact;
    }

    /* somewhat fast case */
    if (xn <= 2 && yn <= 2 && shift <= 2 * FLINT_BITS)
    {
        mp_limb_t t[5], xtmp[4], ytmp[4], yhi, ylo;
        slong fix2;

        xtmp[0] = 0;
        xtmp[1] = 0;

        if (xn == 1)
        {
            xtmp[2] = 0;
            xtmp[3] = xp[0];
        }
        else
        {
            xtmp[2] = xp[0];
            xtmp[3] = xp[1];
        }

        ytmp[0] = 0;
        ytmp[1] = 0;
        ytmp[2] = 0;
        ytmp[3] = 0;

        if (yn == 1)
        {
            yhi = yp[0];
            ylo = 0;
        }
        else
        {
            yhi = yp[1];
            ylo = yp[0];
        }

        shift_limbs = shift / FLINT_BITS;
        shift_bits = shift % FLINT_BITS;

        if (shift_bits == 0)
        {
            ytmp[3 - shift_limbs] = yhi;
            ytmp[2 - shift_limbs] = ylo;
        }
        else
        {
            ytmp[3 - shift_limbs] = yhi >> shift_bits;
            ytmp[2 - shift_limbs] = (yhi << (FLINT_BITS - shift_bits)) | (ylo >> (shift_bits));
            ytmp[1 - shift_limbs] = ylo << (FLINT_BITS - shift_bits);
        }

        if (xsgnbit == ysgnbit)
        {
            t[4] = cy = mpn_add_n(t, xtmp, ytmp, 4);
            fix2 = cy * FLINT_BITS;
            zn = 4 + cy;

            while (t[zn - 1] == 0)
            {
                zn--;
                fix2 -= FLINT_BITS;
            }
        }
        else
        {
#if 1
            cy = mpn_sub_n(t, xtmp, ytmp, 4);

            if (cy)
            {
                mpn_neg(t, t, 4);
                xsgnbit = ysgnbit;
            }

            zn = 4;
            fix2 = 0;

            while (t[zn - 1] == 0)
            {
                zn--;
                fix2 -= FLINT_BITS;

                if (zn == 0)
                {
                    arf_zero(z);
                    return 0;
                }
            }

#else
            int cmp;

            cmp = mpn_cmp(xtmp, ytmp, 4);

            if (cmp > 0)
            {
                mpn_sub_n(t, xtmp, ytmp, 4);
            }
            else if (cmp < 0)
            {
                mpn_sub_n(t, ytmp, xtmp, 4);
                xsgnbit = ysgnbit;
            }
            else
            {
                arf_zero(z);
                return 0;
            }

            fix2 = 0;
            zn = 4;

            while (t[zn - 1] == 0)
            {
                zn--;
                fix2 -= FLINT_BITS;
            }
#endif
        }

        inexact = _arf_set_round_mpn(z, &fix, t, zn, xsgnbit, prec, rnd);
        _fmpz_add_fast(ARF_EXPREF(z), xexp, fix + fix2);
        return inexact;
    }

    /* x +/- eps */
    if (shift > prec + FLINT_BITS + 1 &&
        prec != ARF_PREC_EXACT &&
        shift > (xn + 1) * FLINT_BITS)
    {
        zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        zn_original = zn = FLINT_MAX(zn, xn) + 2;
        shift_limbs = zn - xn;
        alloc = zn + 1;

        ARF_ADD_TMP_ALLOC(tmp, alloc)

        flint_mpn_zero(tmp, shift_limbs);
        flint_mpn_copyi(tmp + shift_limbs, xp, xn);

        if (xsgnbit == ysgnbit)
        {
            tmp[0] = 1;
        }
        else
        {
            mpn_sub_1(tmp, tmp, zn, 1);
            while (tmp[zn-1] == 0)
                zn--;
        }

        _arf_set_round_mpn(z, &fix, tmp, zn, xsgnbit, prec, rnd);
        fix += (zn - zn_original) * FLINT_BITS;
        _fmpz_add_fast(ARF_EXPREF(z), xexp, fix);

        ARF_ADD_TMP_FREE(tmp, alloc)
        return 1;
    }

    shift_limbs = shift / FLINT_BITS;
    shift_bits = shift % FLINT_BITS;

    /* TODO: shift == 0 optimization! (common case) */

    /* w = y shifted */
    wn = yn + (shift_bits != 0);

    /* Size of sum (possibly one more limb than this for carry). */
    zn_original = zn = FLINT_MAX(xn, wn + shift_limbs);
    alloc = zn + 1;

    ARF_ADD_TMP_ALLOC(tmp, alloc)

    /* Low limbs of both operands */
    xbase = zn - xn;
    wbase = zn - shift_limbs - wn;

    /* Shift y to its correct offset in the sum. */
    if (shift_bits == 0)
        flint_mpn_copyi(tmp + wbase, yp, yn);
    else
        tmp[wbase] = cy = mpn_rshift(tmp + wbase + 1, yp, yn, shift_bits);

    if (xsgnbit == ysgnbit)
    {
        if (shift_limbs >= xn)
        {
            /*   XXXXXX
                         WWWWWWW    */
            flint_mpn_zero(tmp + wn, zn - xn - wn);
            flint_mpn_copyi(tmp + zn - xn, xp, xn);
            cy = 0;
        }
        else if (shift_limbs == 0)
        {
            /*  XXXXXX        or    XXXXX
                WWWW                WWWWWWW     */
            if (wn < xn)
            {
                flint_mpn_copyi(tmp, xp, xn - wn);
                cy = mpn_add_n(tmp + wbase, xp + wbase, tmp + wbase, wn);
            }
            else
            {
                cy = mpn_add_n(tmp + xbase, xp, tmp + xbase, xn);
            }
        }
        else if (xbase >= wbase)
        {
            /* XXXXXX         or  XXXXXX
                 WWWWWWW            WWWW    */
            cy = mpn_add_n(tmp + xbase, tmp + xbase, xp, wn - xbase);
            cy = mpn_add_1(tmp + wn, xp + wn - xbase, zn - wn, cy);
        }
        else
        {
            /* XXXXXXXX           todo: merge with above?
                  WWW    */
            flint_mpn_copyi(tmp, xp, wbase);
            cy = mpn_add_n(tmp + wbase, tmp + wbase, xp + wbase, wn);
            cy = mpn_add_1(tmp + zn - shift_limbs, xp + xn - shift_limbs, shift_limbs, cy);
        }

        /* There could be a carry. */
        tmp[zn] = cy;
        zn += cy;
    }
    else
    {
        if (shift_limbs >= xn)
        {
            /*   XXXXXX
                         WWWWWWW    */
            mpn_neg(tmp, tmp, wn);
            flint_mpn_store(tmp + wn, zn - xn - wn, -LIMB_ONE);
            mpn_sub_1(tmp + zn - xn, xp, xn, 1);
        }
        else if (shift_limbs == 0)
        {
            /*  XXXXXX (1)    or    XXXXX   (2)
                WWWW                WWWWWWW     */
            if (wn <= xn) /* (1) */
            {
                if (mpn_cmp(xp + wbase, tmp + wbase, wn) >= 0)
                {
                    flint_mpn_copyi(tmp, xp, wbase);
                    mpn_sub_n(tmp + wbase, xp + wbase, tmp + wbase, wn);
                }
                else
                {
                    xsgnbit = ysgnbit;
                    mpn_sub_n(tmp + wbase, tmp + wbase, xp + wbase, wn);
                    if (wbase != 0)
                    {
                        cy = mpn_neg(tmp, xp, wbase);
                        mpn_sub_1(tmp + wbase, tmp + wbase, wn, cy);
                    }
                }
            }
            else /* (2) */
            {
                if (mpn_cmp(tmp + xbase, xp, xn) >= 0)
                {
                    xsgnbit = ysgnbit;
                    mpn_sub_n(tmp + xbase, tmp + xbase, xp, xn);
                }
                else
                {
                    cy = mpn_neg(tmp, tmp, xbase);
                    mpn_sub_n(tmp + xbase, xp, tmp + xbase, xn);
                    mpn_sub_1(tmp + xbase, tmp + xbase, xn, cy);
                }
            }
        }
        else if (xbase >= wbase)
        {
            /* XXXXXX         or  XXXXXX
                 WWWWWWW            WWWW    */
            if (xbase > 0)
            {
                cy = mpn_sub_n(tmp + xbase, xp, tmp + xbase, wn - xbase);
                cy = mpn_sub_1(tmp + wn, xp + xn - shift_limbs, shift_limbs, cy);
                cy = mpn_neg(tmp, tmp, xbase);
                mpn_sub_1(tmp + xbase, tmp + xbase, xn, cy);
            }
            else
            {
                cy = mpn_sub_n(tmp, xp, tmp, wn);
                cy = mpn_sub_1(tmp + wn, xp + wn, xn - wn, cy);
            }
        }
        else
        {
            /* XXXXXXXX
                  WWW    */
            flint_mpn_copyi(tmp, xp, wbase);
            cy = mpn_sub_n(tmp + wbase, xp + wbase, tmp + wbase, wn);
            mpn_sub_1(tmp + zn - shift_limbs, xp + xn - shift_limbs, shift_limbs, cy);
        }

        /* There could be cancellation. */
        while (zn > 0 && tmp[zn - 1] == 0)
            zn--;
    }

    if (zn == 0)
    {
        arf_zero(z);
        inexact = 0;
    }
    else
    {
        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, xsgnbit, prec, rnd);
        fix += (zn - zn_original) * FLINT_BITS;
        _fmpz_add_fast(ARF_EXPREF(z), xexp, fix);
    }

    ARF_ADD_TMP_FREE(tmp, alloc)
    return inexact;
}

