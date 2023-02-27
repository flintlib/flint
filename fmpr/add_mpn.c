/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

#define ADD_STACK_ALLOC 40
#define ADD_TLS_ALLOC 1000

TLS_PREFIX mp_ptr __add_tmp = NULL;
TLS_PREFIX slong __add_alloc = 0;

void _add_tmp_cleanup(void)
{
    flint_free(__add_tmp);
    __add_tmp = NULL;
    __add_alloc = 0;
}

#define ADD_TMP_ALLOC \
    if (alloc <= ADD_STACK_ALLOC) \
    { \
        tmp = tmp_stack; \
    } \
    else if (alloc <= ADD_TLS_ALLOC) \
    { \
        if (__add_alloc < alloc) \
        { \
            if (__add_alloc == 0) \
            { \
                flint_register_cleanup_function(_add_tmp_cleanup); \
            } \
            __add_tmp = flint_realloc(__add_tmp, sizeof(mp_limb_t) * alloc); \
            __add_alloc = alloc; \
        } \
        tmp = __add_tmp; \
    } \
    else \
    { \
        tmp = flint_malloc(sizeof(mp_limb_t) * alloc); \
    }

#define ADD_TMP_FREE \
    if (alloc > ADD_TLS_ALLOC) \
        flint_free(tmp);

/* computes x + y * 2^shift (optionally negated) */
slong
_fmpr_add_mpn(fmpr_t z,
        mp_srcptr xman, mp_size_t xn, int xsign, const fmpz_t xexp,
        mp_srcptr yman, mp_size_t yn, int ysign, const fmpz_t yexp,
        slong shift, slong prec, fmpr_rnd_t rnd)
{
    slong tn, zn, alloc, ret, shift_bits, shift_limbs;
    int negative;
    mp_limb_t tmp_stack[ADD_STACK_ALLOC];
    mp_limb_t cy;
    mp_ptr tmp, tmp2;

    shift_limbs = shift / FLINT_BITS;
    shift_bits = shift % FLINT_BITS;

    /* x does not overlap with y or the result -- outcome is
       equivalent to adding/subtracting a small number to/from y
       and rounding */
    if (shift > xn * FLINT_BITS &&
        prec != FMPR_PREC_EXACT &&
        xn * FLINT_BITS + prec - (FLINT_BITS * (yn - 1)) < shift)
    {
        zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        zn = FLINT_MAX(zn, yn) + 2;
        shift_limbs = zn - yn;
        alloc = zn;

        ADD_TMP_ALLOC

        flint_mpn_zero(tmp, shift_limbs);
        flint_mpn_copyi(tmp + shift_limbs, yman, yn);

        if (xsign == ysign)
        {
            tmp[0] = 1;
        }
        else
        {
            mpn_sub_1(tmp, tmp, zn, 1);
            while (tmp[zn-1] == 0)
                zn--;
        }

        ret = _fmpr_set_round_mpn(&shift, fmpr_manref(z), tmp, zn, ysign, prec, rnd);
        shift -= shift_limbs * FLINT_BITS;
        fmpz_add_si_inline(fmpr_expref(z), yexp, shift);

        ADD_TMP_FREE

        return ret;
    }

    if (xsign == ysign)
    {
        negative = xsign;

        /* certainly no overlap between x and y * 2^shift

                        XXXXX
               YYYYY         
                    |-------|
                   shift_limbs

          todo: handle huge case efficiently!
        */
        if (shift_limbs >= xn)
        {
            alloc = shift_limbs + yn + 1;

            ADD_TMP_ALLOC

            flint_mpn_copyi(tmp, xman, xn);
            flint_mpn_zero(tmp + xn, shift_limbs - xn);

            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp + shift_limbs, yman, yn);
                zn = shift_limbs + yn;
            }
            else
            {
                cy = mpn_lshift(tmp + shift_limbs, yman, yn, shift_bits);
                tmp[shift_limbs + yn] = cy;
                zn = shift_limbs + yn + (cy != 0);
            }
        }
        else /* may have overlap between x and y * 2^shift */
        {
            slong alloc1, alloc2;

            alloc2 = yn + 1;   /* shifted value of y */
            alloc1 = FLINT_MAX(alloc2 + shift_limbs, xn) + 1; /* space for sum */
            alloc = alloc1 + alloc2;

            ADD_TMP_ALLOC

            tmp2 = tmp + alloc1;

            flint_mpn_copyi(tmp, xman, shift_limbs);

            /* {t, tn} = shift of y (ignoring shift_limbs zero limbs) */
            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp2, yman, yn);
                tn = yn;
            }
            else
            {
                cy = mpn_lshift(tmp2, yman, yn, shift_bits);
                tmp2[yn] = cy;
                tn = yn + (cy != 0);
            }

            if (tn + shift_limbs == xn)
            {
                /*    XXXXXXXX
                      TTTTT        */
                cy = mpn_add_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                tmp[xn] = cy;
                zn = xn + (cy != 0);
            }
            else if (tn + shift_limbs < xn)
            {
                /*      XXXXXXXX
                          TTTTT    */
                cy = mpn_add_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                cy = mpn_add_1(tmp + shift_limbs + tn, xman + shift_limbs + tn, xn - tn - shift_limbs, cy);
                tmp[xn] = cy;
                zn = xn + (cy != 0);
            }
            else
            {
                /*    XXXXXXXX
                   TTTTTTT         */
                cy = mpn_add_n(tmp + shift_limbs, xman + shift_limbs, tmp2, xn - shift_limbs);
                cy = mpn_add_1(tmp + xn, tmp2 + xn - shift_limbs, tn + shift_limbs - xn, cy);
                tmp[tn + shift_limbs] = cy;
                zn = tn + shift_limbs + (cy != 0);
            }
        }
    }
    else
    {
        /* certainly no overlap between x and y * 2^shift

                        XXXXX
               YYYYY         
                    |-------|
                   shift_limbs

          todo: handle huge case efficiently!
        */
        if (shift_limbs >= xn)
        {
            alloc = shift_limbs + yn + 1;

            ADD_TMP_ALLOC

            mpn_neg_n(tmp, xman, xn);

            /* since x is nonzero, there is certainly carry (borrow) */
            flint_mpn_store(tmp + xn, shift_limbs - xn, (mp_limb_t) -1);

            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp + shift_limbs, yman, yn);
                zn = shift_limbs + yn;
            }
            else
            {
                cy = mpn_lshift(tmp + shift_limbs, yman, yn, shift_bits);
                tmp[shift_limbs + yn] = cy;
                zn = shift_limbs + yn + (cy != 0);
            }

            mpn_sub_1(tmp + shift_limbs, tmp + shift_limbs, zn - shift_limbs, 1);

            /* y has larger absolute value, and determines the sign */
            negative = ysign;
        }
        else /* may have overlap between x and y * 2^shift */
        {
            slong alloc1, alloc2;

            alloc2 = yn + 1;   /* shifted value of y */
            alloc1 = FLINT_MAX(alloc2 + shift_limbs, xn); /* space for difference */
            alloc = alloc1 + alloc2;

            ADD_TMP_ALLOC

            tmp2 = tmp + alloc1;

            /* {t, tn} = shift of y (ignoring shift_limbs zero limbs) */
            if (shift_bits == 0)
            {
                flint_mpn_copyi(tmp2, yman, yn);
                tn = yn;
            }
            else
            {
                cy = mpn_lshift(tmp2, yman, yn, shift_bits);
                tmp2[yn] = cy;
                tn = yn + (cy != 0);
            }

            if (tn + shift_limbs == xn)
            {
                int comp;
                /*    XXXXXXXX
                      TTTTT        */

                comp = mpn_cmp(xman + shift_limbs, tmp2, tn);

                if (comp == 0)
                {
                    slong i;

                    for (i = 0; i < shift_limbs; i++)
                    {
                        if (xman[i] != 0)
                        {
                            comp = 1;
                            break;
                        }
                    }

                    /* the result is exactly zero */
                    if (comp == 0)
                    {
                        fmpr_zero(z);
                        ret = FMPR_RESULT_EXACT;
                        goto cleanup;
                    }
                }

                /* x is larger */
                if (comp > 0)
                {
                    flint_mpn_copyi(tmp, xman, shift_limbs);
                    mpn_sub_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                    negative = xsign;
                }
                else
                {
                    mpn_sub_n(tmp + shift_limbs, tmp2, xman + shift_limbs, tn);
                    if (shift_limbs != 0)
                    {
                        cy = mpn_neg_n(tmp, xman, shift_limbs);
                        mpn_sub_1(tmp + shift_limbs, tmp + shift_limbs, tn, cy);
                    }
                    negative = ysign;
                }

                zn = xn;
            }
            else if (tn + shift_limbs < xn)
            {
                /*      XXXXXXXX
                          TTTTT    */

                flint_mpn_copyi(tmp, xman, shift_limbs);
                cy = mpn_sub_n(tmp + shift_limbs, xman + shift_limbs, tmp2, tn);
                mpn_sub_1(tmp + shift_limbs + tn, xman + shift_limbs + tn, xn - tn - shift_limbs, cy);
                zn = xn;
                negative = xsign;
            }
            else
            {
                /*    XXXXXXXX
                   TTTTTTT         */
                cy = mpn_sub_n(tmp + shift_limbs, tmp2, xman + shift_limbs, xn - shift_limbs);
                mpn_sub_1(tmp + xn, tmp2 + xn - shift_limbs, tn + shift_limbs - xn, cy);

                if (shift_limbs != 0)
                {
                    cy = mpn_neg_n(tmp, xman, shift_limbs);
                    mpn_sub_1(tmp + shift_limbs, tmp + shift_limbs, tn, cy);
                }

                zn = shift_limbs + tn;
                negative = ysign;
            }
        }

        /* there might have been cancellation */
        while (tmp[zn-1] == 0)
            zn--;
    }

    ret = _fmpr_set_round_mpn(&shift, fmpr_manref(z), tmp, zn, negative, prec, rnd);

    fmpz_add_si_inline(fmpr_expref(z), xexp, shift);

cleanup:
    ADD_TMP_FREE

    return ret;
}

