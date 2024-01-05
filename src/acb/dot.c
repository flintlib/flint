/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "mpn_extras.h"

/* We need uint64_t instead of mp_limb_t on 32-bit systems for
   safe summation of 30-bit error bounds. */
#include <stdint.h>

/* Add ((a * b) / 2^MAG_BITS) * 2^exp into srad*2^srad_exp.
   Assumes that srad_exp >= exp and that overflow cannot occur.  */
#define RAD_ADDMUL(srad, srad_exp, a, b, exp) \
    do { \
        uint64_t __a, __b; \
        slong __shift; \
        __a = (a); \
        __b = (b); \
        __shift = (srad_exp) - (exp); \
        if (__shift < MAG_BITS) \
            (srad) += (((__a) * (__b)) >> (MAG_BITS + __shift)) + 1; \
        else \
            (srad) += 1; \
    } while (0)

void mag_set_ui_2exp_small(mag_t z, ulong x, slong e);

static void
add_errors(mag_t rad, uint64_t Aerr, slong Aexp, uint64_t Berr, slong Bexp, uint64_t Cerr, slong Cexp)
{
    slong shift;

    if (Aerr && Berr)
    {
        if (Aexp >= Bexp)
        {
            shift = Aexp - Bexp;
            if (shift < 64)
                Aerr = Aerr + (Berr >> shift) + 1;
            else
                Aerr = Aerr + 1;
        }
        else
        {
            shift = Bexp - Aexp;
            if (shift < 64)
                Aerr = Berr + (Aerr >> shift) + 1;
            else
                Aerr = Berr + 1;
            Aexp = Bexp;
        }
    }
    else if (Berr)
    {
        Aerr = Berr;
        Aexp = Bexp;
    }

    if (Aerr && Cerr)
    {
        if (Aexp >= Cexp)
        {
            shift = Aexp - Cexp;
            if (shift < 64)
                Aerr = Aerr + (Cerr >> shift) + 1;
            else
                Aerr = Aerr + 1;
        }
        else
        {
            shift = Cexp - Aexp;
            if (shift < 64)
                Aerr = Cerr + (Aerr >> shift) + 1;
            else
                Aerr = Cerr + 1;
            Aexp = Cexp;
        }
    }
    else if (Cerr)
    {
        Aerr = Cerr;
        Aexp = Cexp;
    }

#if FLINT_BITS == 64
    mag_set_ui_2exp_small(rad, Aerr, Aexp - MAG_BITS);
#else
    mag_set_d(rad, Aerr * (1.0 + 1e-14));
    mag_mul_2exp_si(rad, rad, Aexp - MAG_BITS);
#endif
}

void
_arb_dot_addmul_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn, mp_srcptr yptr, mp_size_t yn,
    int negative, flint_bitcnt_t shift);

void
_arb_dot_add_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn,
    int negative, flint_bitcnt_t shift);

static void
_arb_dot_output(arb_t res, mp_ptr sum, mp_size_t sn, int negative,
    uint64_t serr, slong sum_exp, uint64_t srad, slong srad_exp, slong prec)
{
    slong exp_fix;
    int inexact;

    if (sum[sn - 1] >= LIMB_TOP)
    {
        mpn_neg(sum, sum, sn);
        negative ^= 1;
    }

    exp_fix = 0;

    if (sum[sn - 1] == 0)
    {
        slong sum_exp2;
        mp_size_t sn2;

        sn2 = sn;
        sum_exp2 = sum_exp;

        while (sn2 > 0 && sum[sn2 - 1] == 0)
        {
            sum_exp2 -= FLINT_BITS;
            sn2--;
        }

        if (sn2 == 0)
        {
            arf_zero(arb_midref(res));
            inexact = 0;
        }
        else
        {
            inexact = _arf_set_round_mpn(arb_midref(res), &exp_fix, sum, sn2, negative, prec, ARF_RND_DOWN);
            _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp_fix + sum_exp2);
        }
    }
    else
    {
        if (sn == 2)  /* unnecessary? */
            inexact = _arf_set_round_uiui(arb_midref(res), &exp_fix, sum[1], sum[0], negative, prec, ARF_RND_DOWN);
        else
            inexact = _arf_set_round_mpn(arb_midref(res), &exp_fix, sum, sn, negative, prec, ARF_RND_DOWN);

        _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp_fix + sum_exp);
    }

    add_errors(arb_radref(res),
        inexact << MAG_BITS,
        exp_fix + sum_exp - prec,
        ((uint64_t) serr) << MAG_BITS,
        sum_exp - sn * FLINT_BITS,
        srad,
        srad_exp);
}

/* xxx: don't use surrounding variables */
#define ARB_DOT_ADD(s_sum, s_serr, s_sn, s_sum_exp, s_subtract, xm) \
    if (!arf_is_special(xm)) \
    { \
        mp_srcptr xptr; \
        xexp = ARF_EXP(xm); \
        xn = ARF_SIZE(xm); \
        xnegative = ARF_SGNBIT(xm); \
        shift = s_sum_exp - xexp; \
        if (shift >= s_sn * FLINT_BITS) \
        { \
            s_serr++; \
        } \
        else \
        { \
            xptr = (xn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(xm) : ARF_PTR_D(xm); \
            _arb_dot_add_generic(s_sum, &s_serr, tmp, s_sn, xptr, xn, xnegative ^ s_subtract, shift); \
        } \
    } \

/* xxx: don't use surrounding variables */
#define ARB_DOT_ADD_RAD(s_srad, s_srad_exp, xr) \
    if (!mag_is_special(xr)) \
    { \
        xrad = MAG_MAN(xr); \
        xrexp = MAG_EXP(xr); \
        shift = s_srad_exp - xrexp; \
        if (shift < 64) \
            s_srad += (xrad >> shift) + 1; \
        else \
            s_srad++; \
    }

static void
_arf_complex_mul_gauss(arf_t e, arf_t f, const arf_t a, const arf_t b,
                                         const arf_t c, const arf_t d)
{
    mp_srcptr ap, bp, cp, dp;
    int asgn, bsgn, csgn, dsgn;
    mp_size_t an, bn, cn, dn;
    slong aexp, bexp, cexp, dexp;
    fmpz texp, uexp;

    fmpz_t za, zb, zc, zd, t, u, v;
    slong abot, bbot, cbot, dbot;

    ARF_GET_MPN_READONLY(ap, an, a);
    asgn = ARF_SGNBIT(a);
    aexp = ARF_EXP(a);

    ARF_GET_MPN_READONLY(bp, bn, b);
    bsgn = ARF_SGNBIT(b);
    bexp = ARF_EXP(b);

    ARF_GET_MPN_READONLY(cp, cn, c);
    csgn = ARF_SGNBIT(c);
    cexp = ARF_EXP(c);

    ARF_GET_MPN_READONLY(dp, dn, d);
    dsgn = ARF_SGNBIT(d);
    dexp = ARF_EXP(d);

    /* Gauss multiplication
        e = ac - bd
        f = (a+b)(c+d) - ac - bd */

    abot = aexp - an * FLINT_BITS;
    bbot = bexp - bn * FLINT_BITS;
    cbot = cexp - cn * FLINT_BITS;
    dbot = dexp - dn * FLINT_BITS;

    texp = FLINT_MIN(abot, bbot);
    uexp = FLINT_MIN(cbot, dbot);

    fmpz_init(za);
    fmpz_init(zb);
    fmpz_init(zc);
    fmpz_init(zd);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(v);

    fmpz_lshift_mpn(za, ap, an, asgn, abot - texp);
    fmpz_lshift_mpn(zb, bp, bn, bsgn, bbot - texp);
    fmpz_lshift_mpn(zc, cp, cn, csgn, cbot - uexp);
    fmpz_lshift_mpn(zd, dp, dn, dsgn, dbot - uexp);

    fmpz_add(t, za, zb);
    fmpz_add(v, zc, zd);
    fmpz_mul(u, t, v);
    fmpz_mul(t, za, zc);
    fmpz_mul(v, zb, zd);
    fmpz_sub(u, u, t);
    fmpz_sub(u, u, v);
    fmpz_sub(t, t, v);

    texp += uexp;
    arf_set_fmpz_2exp(e, t, &texp);
    arf_set_fmpz_2exp(f, u, &texp);

    fmpz_clear(za);
    fmpz_clear(zb);
    fmpz_clear(zc);
    fmpz_clear(zd);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(v);
}

/* TODO: this could be much lower, but it's currently competing
   against mulhigh in the Karatsuba range. */
FLINT_DLL slong acb_dot_gauss_dot_cutoff = 128;
#define GAUSS_CUTOFF acb_dot_gauss_dot_cutoff

void
acb_dot(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)
{
    slong i, j, padding, extend;
    slong xexp, yexp, exp;
    slong re_nonzero, im_nonzero;
    slong re_max_exp, re_min_exp, re_sum_exp;
    slong im_max_exp, im_min_exp, im_sum_exp;
    slong re_srad_exp, re_max_rad_exp;
    slong im_srad_exp, im_max_rad_exp;
    slong re_prec, im_prec;
    slong xrexp, yrexp;
    int xnegative, ynegative;
    mp_size_t xn, yn, re_sn, im_sn, alloc;
    flint_bitcnt_t shift;
    arb_srcptr xi, yi;
    arf_srcptr xm, ym;
    mag_srcptr xr, yr;
    mp_limb_t xtop, ytop;
    mp_limb_t xrad, yrad;
    mp_limb_t re_serr, im_serr;   /* Sum over arithmetic errors */
    uint64_t re_srad, im_srad;    /* Sum over propagated errors */
    mp_ptr tmp, re_sum, im_sum;   /* Workspace */
    slong xoff, yoff;
    char * use_gauss;
    ARF_ADD_TMP_DECL;

    /* todo: fast fma and fmma (len=2) code */
    if (len <= 1)
    {
        if (initial == NULL)
        {
            if (len <= 0)
                acb_zero(res);
            else
            {
                acb_mul(res, x, y, prec);
                if (subtract)
                    acb_neg(res, res);
            }
            return;
        }
        else if (len <= 0)
        {
            acb_set_round(res, initial, prec);
            return;
        }
    }

    /* Number of nonzero midpoint terms in sum. */
    re_nonzero = 0;
    im_nonzero = 0;

    /* Terms are bounded by 2^max_exp (with WORD_MIN = -infty) */
    re_max_exp = WORD_MIN;
    im_max_exp = WORD_MIN;

    /* Propagated error terms are bounded by 2^max_rad_exp */
    re_max_rad_exp = WORD_MIN;
    im_max_rad_exp = WORD_MIN;

    /* Used to reduce the precision. */
    re_min_exp = WORD_MAX;
    im_min_exp = WORD_MAX;

    /* Account for the initial term. */
    if (initial != NULL)
    {
        if (!ARB_IS_LAGOM(acb_realref(initial)) || !ARB_IS_LAGOM(acb_imagref(initial)))
        {
            acb_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec);
            return;
        }

        xm = arb_midref(acb_realref(initial));
        xr = arb_radref(acb_realref(initial));

        if (!arf_is_special(xm))
        {
            re_max_exp = ARF_EXP(xm);
            re_nonzero++;

            if (prec > 2 * FLINT_BITS)
                re_min_exp = ARF_EXP(xm) - ARF_SIZE(xm) * FLINT_BITS;
        }

        if (!mag_is_special(xr))
            re_max_rad_exp = MAG_EXP(xr);

        xm = arb_midref(acb_imagref(initial));
        xr = arb_radref(acb_imagref(initial));

        if (!arf_is_special(xm))
        {
            im_max_exp = ARF_EXP(xm);
            im_nonzero++;

            if (prec > 2 * FLINT_BITS)
                im_min_exp = ARF_EXP(xm) - ARF_SIZE(xm) * FLINT_BITS;
        }

        if (!mag_is_special(xr))
            im_max_rad_exp = MAG_EXP(xr);
    }

    for (xoff = 0; xoff < 2; xoff++)
    {
        for (yoff = 0; yoff < 2; yoff++)
        {
            slong nonzero, max_exp, min_exp, max_rad_exp;

            if (xoff == yoff)
            {
                nonzero = re_nonzero;
                max_exp = re_max_exp;
                min_exp = re_min_exp;
                max_rad_exp = re_max_rad_exp;
            }
            else
            {
                nonzero = im_nonzero;
                max_exp = im_max_exp;
                min_exp = im_min_exp;
                max_rad_exp = im_max_rad_exp;
            }

            /* Determine maximum exponents for the main sum and the radius sum. */
            for (i = 0; i < len; i++)
            {
                xi = ((arb_srcptr) x) + 2 * i * xstep + xoff;
                yi = ((arb_srcptr) y) + 2 * i * ystep + yoff;

                /* Fallback for huge exponents or non-finite values. */
                if (!ARB_IS_LAGOM(xi) || !ARB_IS_LAGOM(yi))
                {
                    acb_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec);
                    return;
                }

                xm = arb_midref(xi);
                ym = arb_midref(yi);
                xr = arb_radref(xi);
                yr = arb_radref(yi);

                /* (xm+xr)(ym+yr) = xm ym + [xr ym + xm yr + xr yr] */
                if (!arf_is_special(xm))
                {
                    xexp = ARF_EXP(xm);

                    if (!arf_is_special(ym))
                    {
                        yexp = ARF_EXP(ym);

                        max_exp = FLINT_MAX(max_exp, xexp + yexp);
                        nonzero++;

                        if (prec > 2 * FLINT_BITS)
                        {
                            slong bot;
                            bot = (xexp + yexp) - (ARF_SIZE(xm) + ARF_SIZE(ym)) * FLINT_BITS;
                            min_exp = FLINT_MIN(min_exp, bot);
                        }

                        if (!mag_is_special(xr))
                        {
                            xrexp = MAG_EXP(xr);
                            max_rad_exp = FLINT_MAX(max_rad_exp, yexp + xrexp);

                            if (!mag_is_special(yr))
                            {
                                yrexp = MAG_EXP(yr);
                                max_rad_exp = FLINT_MAX(max_rad_exp, xexp + yrexp);
                                max_rad_exp = FLINT_MAX(max_rad_exp, xrexp + yrexp);
                            }
                        }
                        else
                        {
                            if (!mag_is_special(yr))
                            {
                                yrexp = MAG_EXP(yr);
                                max_rad_exp = FLINT_MAX(max_rad_exp, xexp + yrexp);
                            }
                        }
                    }
                    else  /* if y = 0, something can happen only if yr != 0 */
                    {
                        if (!mag_is_special(yr))
                        {
                            yrexp = MAG_EXP(yr);
                            max_rad_exp = FLINT_MAX(max_rad_exp, xexp + yrexp);

                            if (!mag_is_special(xr))
                            {
                                xrexp = MAG_EXP(xr);
                                max_rad_exp = FLINT_MAX(max_rad_exp, xrexp + yrexp);
                            }
                        }
                    }
                }
                else  /* if x = 0, something can happen only if xr != 0 */
                {
                    if (!mag_is_special(xr))
                    {
                        xrexp = MAG_EXP(xr);

                        if (!arf_is_special(ym))
                        {
                            yexp = ARF_EXP(ym);
                            max_rad_exp = FLINT_MAX(max_rad_exp, xrexp + yexp);
                        }

                        if (!mag_is_special(yr))
                        {
                            yrexp = MAG_EXP(yr);
                            max_rad_exp = FLINT_MAX(max_rad_exp, xrexp + yrexp);
                        }
                    }
                }
            }

            if (xoff == yoff)
            {
                re_nonzero = nonzero;
                re_max_exp = max_exp;
                re_min_exp = min_exp;
                re_max_rad_exp = max_rad_exp;
            }
            else
            {
                im_nonzero = nonzero;
                im_max_exp = max_exp;
                im_min_exp = min_exp;
                im_max_rad_exp = max_rad_exp;
            }
        }
    }

    re_prec = prec;
    im_prec = prec;

    if (re_max_exp == WORD_MIN && re_max_rad_exp == WORD_MIN &&
        im_max_exp == WORD_MIN && im_max_rad_exp == WORD_MIN)
    {
        acb_zero(res);
        return;
    }

    /* The midpoint sum is zero. */
    if (re_max_exp == WORD_MIN)
    {
        re_prec = 2;
    }
    else
    {
        if (re_max_rad_exp != WORD_MIN)
            re_prec = FLINT_MIN(re_prec, re_max_exp - re_max_rad_exp + MAG_BITS);
        if (re_min_exp != WORD_MAX)
            re_prec = FLINT_MIN(re_prec, re_max_exp - re_min_exp + MAG_BITS);
        re_prec = FLINT_MAX(re_prec, 2);
    }

    if (im_max_exp == WORD_MIN)
    {
        im_prec = 2;
    }
    else
    {
        if (im_max_rad_exp != WORD_MIN)
            im_prec = FLINT_MIN(im_prec, im_max_exp - im_max_rad_exp + MAG_BITS);
        if (re_min_exp != WORD_MAX)
            im_prec = FLINT_MIN(im_prec, im_max_exp - im_min_exp + MAG_BITS);
        im_prec = FLINT_MAX(im_prec, 2);
    }

    extend = FLINT_BIT_COUNT(re_nonzero) + 1;
    padding = 4 + FLINT_BIT_COUNT(len);
    re_sn = (re_prec + extend + padding + FLINT_BITS - 1) / FLINT_BITS;
    re_sn = FLINT_MAX(re_sn, 2);
    re_sum_exp = re_max_exp + extend;

    extend = FLINT_BIT_COUNT(im_nonzero) + 1;
    padding = 4 + FLINT_BIT_COUNT(len);
    im_sn = (im_prec + extend + padding + FLINT_BITS - 1) / FLINT_BITS;
    im_sn = FLINT_MAX(im_sn, 2);
    im_sum_exp = im_max_exp + extend;

    /* We need sn + 1 limb for the sum (sn limbs + 1 dummy limb
       for carry or borrow that avoids an extra branch). We need
       2 * (sn + 2) limbs to store the product of two numbers
       with up to (sn + 2) limbs, plus 1 extra limb for shifting
       the product. */
    alloc = (re_sn + 1) + (im_sn + 1) + 2 * (FLINT_MAX(re_sn, im_sn) + 2) + 1;
    ARF_ADD_TMP_ALLOC(re_sum, alloc)
    im_sum = re_sum + (re_sn + 1);
    tmp = im_sum + (im_sn + 1);

    /* Sum of propagated errors. */
    re_srad_exp = re_max_rad_exp;
    re_srad = 0;
    im_srad_exp = im_max_rad_exp;
    im_srad = 0;

    /* Set sum to 0 */
    re_serr = 0;
    for (j = 0; j < re_sn + 1; j++)
        re_sum[j] = 0;
    im_serr = 0;
    for (j = 0; j < im_sn + 1; j++)
        im_sum[j] = 0;

    if (initial != NULL)
    {
        xm = arb_midref(acb_realref(initial));
        xr = arb_radref(acb_realref(initial));

        ARB_DOT_ADD(re_sum, re_serr, re_sn, re_sum_exp, subtract, xm);
        ARB_DOT_ADD_RAD(re_srad, re_srad_exp, xr);

        xm = arb_midref(acb_imagref(initial));
        xr = arb_radref(acb_imagref(initial));

        ARB_DOT_ADD(im_sum, im_serr, im_sn, im_sum_exp, subtract, xm);
        ARB_DOT_ADD_RAD(im_srad, im_srad_exp, xr);
    }

    /*
    Look for terms to process using the Gauss multiplication formula.
    If any such terms are found, we mask the ith entry in use_gauss
    so that they will be skipped in the main loop.
    Important: the cutoffs must be such that the fast case
    (xn <= 2, yn <= 2, sn <= 3) is not hit below and the mask
    check is done.

    The cutoffs below are not optimal in the generic case; also, it
    would be nicer to have both mulhigh and Gauss here. A more elegant
    solution would be to write a fallback version of acb_dot_simple
    where acb_addmul does the right thing.
    */
    use_gauss = NULL;

    if (re_prec >= GAUSS_CUTOFF * FLINT_BITS &&
        im_prec >= GAUSS_CUTOFF * FLINT_BITS)
    {
        arf_t e, f;

        for (i = 0; i < len; i++)
        {
            arb_srcptr ai, bi, ci, di;
            mp_size_t an, bn, cn, dn;
            slong aexp, bexp, cexp, dexp;

            ai = ((arb_srcptr) x) + 2 * i * xstep;
            bi = ((arb_srcptr) x) + 2 * i * xstep + 1;
            ci = ((arb_srcptr) y) + 2 * i * ystep;
            di = ((arb_srcptr) y) + 2 * i * ystep + 1;

            an = ARF_SIZE(arb_midref(ai));
            bn = ARF_SIZE(arb_midref(bi));
            cn = ARF_SIZE(arb_midref(ci));
            dn = ARF_SIZE(arb_midref(di));

            aexp = ARF_EXP(arb_midref(ai));
            bexp = ARF_EXP(arb_midref(bi));
            cexp = ARF_EXP(arb_midref(ci));
            dexp = ARF_EXP(arb_midref(di));

            if (an >= GAUSS_CUTOFF && bn >= GAUSS_CUTOFF &&
                bn >= GAUSS_CUTOFF && cn >= GAUSS_CUTOFF &&
                FLINT_ABS(an - bn) <= 2 &&
                FLINT_ABS(cn - dn) <= 2 &&
                FLINT_ABS(aexp - bexp) <= 64 &&
                FLINT_ABS(cexp - dexp) <= 64 &&
                re_sum_exp - (aexp + cexp) < 0.1 * re_prec &&
                im_sum_exp - (aexp + dexp) < 0.1 * im_prec &&
                an + cn < 2.2 * re_sn && an + dn < 2.2 * im_sn)
            {
                if (use_gauss == NULL)
                {
                    use_gauss = flint_calloc(len, sizeof(char));
                    arf_init(e);
                    arf_init(f);
                }

                use_gauss[i] = 1;
                _arf_complex_mul_gauss(e, f, arb_midref(ai), arb_midref(bi), arb_midref(ci), arb_midref(di));
                ARB_DOT_ADD(re_sum, re_serr, re_sn, re_sum_exp, 0, e);
                ARB_DOT_ADD(im_sum, im_serr, im_sn, im_sum_exp, 0, f);
            }
        }

        if (use_gauss != NULL)
        {
            arf_clear(e);
            arf_clear(f);
        }
    }

    for (xoff = 0; xoff < 2; xoff++)
    {
        for (yoff = 0; yoff < 2; yoff++)
        {
            slong sum_exp, srad_exp;
            mp_ptr sum;
            mp_size_t sn;
            mp_limb_t serr;
            uint64_t srad;
            int flipsign;

            if (xoff == yoff)
            {
                sum_exp = re_sum_exp;
                srad_exp = re_srad_exp;
                sum = re_sum;
                sn = re_sn;
                if (re_max_exp == WORD_MIN && re_max_rad_exp == WORD_MIN)
                    continue;
            }
            else
            {
                sum_exp = im_sum_exp;
                srad_exp = im_srad_exp;
                sum = im_sum;
                sn = im_sn;
                if (im_max_exp == WORD_MIN && im_max_rad_exp == WORD_MIN)
                    continue;
            }

            serr = 0;
            srad = 0;
            flipsign = (xoff + yoff == 2);

            for (i = 0; i < len; i++)
            {
                xi = ((arb_srcptr) x) + 2 * i * xstep + xoff;
                yi = ((arb_srcptr) y) + 2 * i * ystep + yoff;

                xm = arb_midref(xi);
                ym = arb_midref(yi);
                xr = arb_radref(xi);
                yr = arb_radref(yi);

                /* The midpoints of x[i] and y[i] are both nonzero. */
                if (!arf_is_special(xm) && !arf_is_special(ym))
                {
                    xexp = ARF_EXP(xm);
                    xn = ARF_SIZE(xm);
                    xnegative = ARF_SGNBIT(xm);

                    yexp = ARF_EXP(ym);
                    yn = ARF_SIZE(ym);
                    ynegative = ARF_SGNBIT(ym);

                    exp = xexp + yexp;
                    shift = sum_exp - exp;

                    if (shift >= sn * FLINT_BITS)
                    {
                        /* We may yet need the top limbs for bounds. */
                        ARF_GET_TOP_LIMB(xtop, xm);
                        ARF_GET_TOP_LIMB(ytop, ym);
                        serr++;
                    }
                    else if (xn <= 2 && yn <= 2 && sn <= 3)
                    {
                        mp_limb_t x1, x0, y1, y0;
                        mp_limb_t u3, u2, u1, u0;

                        if (xn == 1 && yn == 1)
                        {
                            xtop = ARF_NOPTR_D(xm)[0];
                            ytop = ARF_NOPTR_D(ym)[0];
                            umul_ppmm(u3, u2, xtop, ytop);
                            u1 = u0 = 0;
                        }
                        else if (xn == 2 && yn == 2)
                        {
                            x0 = ARF_NOPTR_D(xm)[0];
                            x1 = ARF_NOPTR_D(xm)[1];
                            y0 = ARF_NOPTR_D(ym)[0];
                            y1 = ARF_NOPTR_D(ym)[1];
                            xtop = x1;
                            ytop = y1;
                            FLINT_MPN_MUL_2X2(u3, u2, u1, u0, x1, x0, y1, y0);
                        }
                        else if (xn == 1)
                        {
                            x0 = ARF_NOPTR_D(xm)[0];
                            y0 = ARF_NOPTR_D(ym)[0];
                            y1 = ARF_NOPTR_D(ym)[1];
                            xtop = x0;
                            ytop = y1;
                            FLINT_MPN_MUL_2X1(u3, u2, u1, y1, y0, x0);
                            u0 = 0;
                        }
                        else
                        {
                            x0 = ARF_NOPTR_D(xm)[0];
                            x1 = ARF_NOPTR_D(xm)[1];
                            y0 = ARF_NOPTR_D(ym)[0];
                            xtop = x1;
                            ytop = y0;
                            FLINT_MPN_MUL_2X1(u3, u2, u1, x1, x0, y0);
                            u0 = 0;
                        }

                        if (sn == 2)
                        {
                            if (shift < FLINT_BITS)
                            {
                                serr += ((u2 << (FLINT_BITS - shift)) != 0) || (u1 != 0) || (u0 != 0);
                                u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                                u3 = (u3 >> shift);
                            }
                            else if (shift == FLINT_BITS)
                            {
                                serr += (u2 != 0) || (u1 != 0) || (u0 != 0);
                                u2 = u3;
                                u3 = 0;
                            }
                            else /* FLINT_BITS < shift < 2 * FLINT_BITS */
                            {
                                serr += ((u3 << (2 * FLINT_BITS - shift)) != 0) || (u2 != 0) || (u1 != 0) || (u0 != 0);
                                u2 = (u3 >> (shift - FLINT_BITS));
                                u3 = 0;
                            }

                            if (xnegative ^ ynegative ^ flipsign)
                                sub_ddmmss(sum[1], sum[0], sum[1], sum[0], u3, u2);
                            else
                                add_ssaaaa(sum[1], sum[0], sum[1], sum[0], u3, u2);
                        }
                        else if (sn == 3)
                        {
                            if (shift < FLINT_BITS)
                            {
                                serr += ((u1 << (FLINT_BITS - shift)) != 0) || (u0 != 0);
                                u1 = (u1 >> shift) | (u2 << (FLINT_BITS - shift));
                                u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                                u3 = (u3 >> shift);
                            }
                            else if (shift == FLINT_BITS)
                            {
                                serr += (u1 != 0) || (u0 != 0);
                                u1 = u2;
                                u2 = u3;
                                u3 = 0;
                            }
                            else if (shift < 2 * FLINT_BITS)
                            {
                                serr += ((u2 << (2 * FLINT_BITS - shift)) != 0) || (u1 != 0) || (u0 != 0);
                                u1 = (u3 << (2 * FLINT_BITS - shift)) | (u2 >> (shift - FLINT_BITS));
                                u2 = (u3 >> (shift - FLINT_BITS));
                                u3 = 0;
                            }
                            else if (shift == 2 * FLINT_BITS)
                            {
                                serr += (u2 != 0) || (u1 != 0) || (u0 != 0);
                                u1 = u3;
                                u2 = 0;
                                u3 = 0;
                            }
                            else  /* 2 * FLINT_BITS < shift < 3 * FLINT_BITS */
                            {
                                serr += ((u3 << (3 * FLINT_BITS - shift)) != 0) || (u2 != 0) || (u1 != 0) || (u0 != 0);
                                u1 = (u3 >> (shift - 2 * FLINT_BITS));
                                u2 = 0;
                                u3 = 0;
                            }

                            if (xnegative ^ ynegative ^ flipsign)
                                sub_dddmmmsss(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                            else
                                add_sssaaaaaa(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                        }
                    }
                    else
                    {
                        mp_srcptr xptr, yptr;

                        xptr = (xn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(xm) : ARF_PTR_D(xm);
                        yptr = (yn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(ym) : ARF_PTR_D(ym);

                        xtop = xptr[xn - 1];
                        ytop = yptr[yn - 1];

                        if (use_gauss == NULL || use_gauss[i] == 0)
                            _arb_dot_addmul_generic(sum, &serr, tmp, sn, xptr, xn, yptr, yn, xnegative ^ ynegative ^ flipsign, shift);
                    }

                    xrad = MAG_MAN(xr);
                    yrad = MAG_MAN(yr);

                    if (xrad != 0 && yrad != 0)
                    {
                        xrexp = MAG_EXP(xr);
                        yrexp = MAG_EXP(yr);

                        RAD_ADDMUL(srad, srad_exp, (xtop >> (FLINT_BITS - MAG_BITS)) + 1, yrad, xexp + yrexp);
                        RAD_ADDMUL(srad, srad_exp, (ytop >> (FLINT_BITS - MAG_BITS)) + 1, xrad, yexp + xrexp);
                        RAD_ADDMUL(srad, srad_exp, xrad, yrad, xrexp + yrexp);
                    }
                    else if (xrad != 0)
                    {
                        xrexp = MAG_EXP(xr);
                        RAD_ADDMUL(srad, srad_exp, (ytop >> (FLINT_BITS - MAG_BITS)) + 1, xrad, yexp + xrexp);
                    }
                    else if (yrad != 0)
                    {
                        yrexp = MAG_EXP(yr);
                        RAD_ADDMUL(srad, srad_exp, (xtop >> (FLINT_BITS - MAG_BITS)) + 1, yrad, xexp + yrexp);
                    }
                }
                else
                {
                    xrad = MAG_MAN(xr);
                    yrad = MAG_MAN(yr);

                    xexp = ARF_EXP(xm);
                    yexp = ARF_EXP(ym);

                    xrexp = MAG_EXP(xr);
                    yrexp = MAG_EXP(yr);

                    /* (xm+xr)(ym+yr) = xm ym + [xm yr + ym xr + xr yr] */
                    if (yrad && !arf_is_special(xm))
                    {
                        ARF_GET_TOP_LIMB(xtop, xm);
                        RAD_ADDMUL(srad, srad_exp, (xtop >> (FLINT_BITS - MAG_BITS)) + 1, yrad, xexp + yrexp);
                    }

                    if (xrad && !arf_is_special(ym))
                    {
                        ARF_GET_TOP_LIMB(ytop, ym);
                        RAD_ADDMUL(srad, srad_exp, (ytop >> (FLINT_BITS - MAG_BITS)) + 1, xrad, yexp + xrexp);
                    }

                    if (xrad && yrad)
                    {
                        RAD_ADDMUL(srad, srad_exp, xrad, yrad, xrexp + yrexp);
                    }
                }
            }

            if (xoff == yoff)
            {
                re_serr += serr;
                re_srad += srad;
            }
            else
            {
                im_serr += serr;
                im_srad += srad;
            }
        }
    }

    _arb_dot_output(acb_realref(res), re_sum, re_sn, subtract, re_serr, re_sum_exp, re_srad, re_srad_exp, re_prec);
    _arb_dot_output(acb_imagref(res), im_sum, im_sn, subtract, im_serr, im_sum_exp, im_srad, im_srad_exp, im_prec);

    ARF_ADD_TMP_FREE(re_sum, alloc);
    if (use_gauss != NULL)
        flint_free(use_gauss);
}
