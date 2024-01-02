/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acf.h"
#include "mpn_extras.h"

/* We need uint64_t instead of mp_limb_t on 32-bit systems for
   safe summation of 30-bit error bounds. */
#include <stdint.h>

void
_arb_dot_addmul_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn, mp_srcptr yptr, mp_size_t yn,
    int negative, flint_bitcnt_t shift);

void
_arb_dot_add_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn,
    int negative, flint_bitcnt_t shift);

static void
_arb_dot_output(arf_t res, mp_ptr sum, mp_size_t sn, int negative,
    slong sum_exp, slong prec, arf_rnd_t rnd)
{
    slong exp_fix;

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
            arf_zero(res);
        }
        else
        {
            _arf_set_round_mpn(res, &exp_fix, sum, sn2, negative, prec, rnd);
            _fmpz_set_si_small(ARF_EXPREF(res), exp_fix + sum_exp2);
        }
    }
    else
    {
        if (sn == 2)  /* unnecessary? */
            _arf_set_round_uiui(res, &exp_fix, sum[1], sum[0], negative, prec, rnd);
        else
            _arf_set_round_mpn(res, &exp_fix, sum, sn, negative, prec, rnd);

        _fmpz_set_si_small(ARF_EXPREF(res), exp_fix + sum_exp);
    }
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
        } \
        else \
        { \
            xptr = (xn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(xm) : ARF_PTR_D(xm); \
            _arb_dot_add_generic(s_sum, &s_serr, tmp, s_sn, xptr, xn, xnegative ^ s_subtract, shift); \
        } \
    } \

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

    /* Karstsuba multiplication
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

FLINT_DLL extern slong acb_dot_gauss_dot_cutoff;
#define GAUSS_CUTOFF acb_dot_gauss_dot_cutoff

void
acf_approx_dot_simple(acf_t res, const acf_t initial, int subtract,
    acf_srcptr x, slong xstep, acf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
        {
            arf_zero(acf_realref(res));
            arf_zero(acf_imagref(res));
        }
        else
        {
            arf_set_round(acf_realref(res), acf_realref(initial), prec, rnd);
            arf_set_round(acf_imagref(res), acf_imagref(initial), prec, rnd);
        }
        return;
    }

    if (initial == NULL && len == 1)
    {
        arf_complex_mul(acf_realref(res),
                        acf_imagref(res),
                        acf_realref(x),
                        acf_imagref(x),
                        acf_realref(y),
                        acf_imagref(y), prec, rnd);
    }
    else
    {
        arf_t e, f;

        arf_init(e);
        arf_init(f);

        if (initial != NULL)
        {
            if (subtract)
            {
                arf_neg(acf_realref(res), acf_realref(initial));
                arf_neg(acf_imagref(res), acf_imagref(initial));
            }
            else
            {
                arf_set(acf_realref(res), acf_realref(initial));
                arf_set(acf_imagref(res), acf_imagref(initial));
            }
        }

        for (i = 0; i < len; i++)
        {
            arf_complex_mul(e, f,
                            (acf_realref(x + i * xstep)),
                            (acf_imagref(x + i * xstep)),
                            (acf_realref(y + i * ystep)),
                            (acf_imagref(y + i * ystep)), prec, rnd);


            if (i == 0 && initial == NULL)
            {
                arf_set(acf_realref(res), e);
                arf_set(acf_imagref(res), f);
            }
            else
            {
                arf_add(acf_realref(res), acf_realref(res), e, prec, rnd);
                arf_add(acf_imagref(res), acf_imagref(res), f, prec, rnd);
            }
        }

        arf_clear(e);
        arf_clear(f);
    }

    if (subtract)
    {
        arf_neg(acf_realref(res), acf_realref(res));
        arf_neg(acf_imagref(res), acf_imagref(res));
    }
}

void
acf_approx_dot(acf_t res, const acf_t initial, int subtract, acf_srcptr x, slong xstep, acf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)
{
    slong i, j, padding, extend;
    slong xexp, yexp, exp;
    slong re_nonzero, im_nonzero;
    slong re_max_exp, re_min_exp, re_sum_exp;
    slong im_max_exp, im_min_exp, im_sum_exp;
    slong re_prec, im_prec;
    int xnegative, ynegative;
    mp_size_t xn, yn, re_sn, im_sn, alloc;
    flint_bitcnt_t shift;
    arf_srcptr xi, yi;
    arf_srcptr xm, ym;
    mp_limb_t re_serr, im_serr;   /* Sum over arithmetic errors */
    mp_ptr tmp, re_sum, im_sum;   /* Workspace */
    slong xoff, yoff;
    char * use_gauss;
    ARF_ADD_TMP_DECL;

    /* todo: fast fma and fmma (len=2) code */
    if (len <= 1)
    {
        acf_approx_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec, rnd);
        return;
    }

    /* Number of nonzero midpoint terms in sum. */
    re_nonzero = 0;
    im_nonzero = 0;

    /* Terms are bounded by 2^max_exp (with WORD_MIN = -infty) */
    re_max_exp = WORD_MIN;
    im_max_exp = WORD_MIN;

    /* Used to reduce the precision. */
    re_min_exp = WORD_MAX;
    im_min_exp = WORD_MAX;

    /* Account for the initial term. */
    if (initial != NULL)
    {
        if (!ARF_IS_LAGOM(acf_realref(initial)) || !ARF_IS_LAGOM(acf_imagref(initial)))
        {
            acf_approx_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec, rnd);
            return;
        }

        xm = acf_realref(initial);

        if (!arf_is_special(xm))
        {
            re_max_exp = ARF_EXP(xm);
            re_nonzero++;

            if (prec > 2 * FLINT_BITS)
                re_min_exp = ARF_EXP(xm) - ARF_SIZE(xm) * FLINT_BITS;
        }

        xm = acf_imagref(initial);

        if (!arf_is_special(xm))
        {
            im_max_exp = ARF_EXP(xm);
            im_nonzero++;

            if (prec > 2 * FLINT_BITS)
                im_min_exp = ARF_EXP(xm) - ARF_SIZE(xm) * FLINT_BITS;
        }
    }

    for (xoff = 0; xoff < 2; xoff++)
    {
        for (yoff = 0; yoff < 2; yoff++)
        {
            slong nonzero, max_exp, min_exp;

            if (xoff == yoff)
            {
                nonzero = re_nonzero;
                max_exp = re_max_exp;
                min_exp = re_min_exp;
            }
            else
            {
                nonzero = im_nonzero;
                max_exp = im_max_exp;
                min_exp = im_min_exp;
            }

            /* Determine maximum exponents for the main sum and the radius sum. */
            for (i = 0; i < len; i++)
            {
                xi = ((arf_srcptr) x) + 2 * i * xstep + xoff;
                yi = ((arf_srcptr) y) + 2 * i * ystep + yoff;

                /* Fallback for huge exponents or non-finite values. */
                if (!ARF_IS_LAGOM(xi) || !ARF_IS_LAGOM(yi))
                {
                    acf_approx_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec, rnd);
                    return;
                }

                xm = xi;
                ym = yi;

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
                    }
                }
            }

            if (xoff == yoff)
            {
                re_nonzero = nonzero;
                re_max_exp = max_exp;
                re_min_exp = min_exp;
            }
            else
            {
                im_nonzero = nonzero;
                im_max_exp = max_exp;
                im_min_exp = min_exp;
            }
        }
    }

    re_prec = prec;
    im_prec = prec;

    if (re_max_exp == WORD_MIN && im_max_exp == WORD_MIN)
    {
        arf_zero(acf_realref(res));
        arf_zero(acf_imagref(res));
        return;
    }

    /* The midpoint sum is zero. */
    if (re_max_exp == WORD_MIN)
    {
        re_prec = 2;
    }
    else
    {
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

    /* Set sum to 0 */
    re_serr = 0;
    for (j = 0; j < re_sn + 1; j++)
        re_sum[j] = 0;
    im_serr = 0;
    for (j = 0; j < im_sn + 1; j++)
        im_sum[j] = 0;

    if (initial != NULL)
    {
        xm = acf_realref(initial);

        ARB_DOT_ADD(re_sum, re_serr, re_sn, re_sum_exp, subtract, xm);

        xm = acf_imagref(initial);

        ARB_DOT_ADD(im_sum, im_serr, im_sn, im_sum_exp, subtract, xm);
    }

    use_gauss = NULL;

    if (re_prec >= GAUSS_CUTOFF * FLINT_BITS &&
        im_prec >= GAUSS_CUTOFF * FLINT_BITS)
    {
        arf_t e, f;

        for (i = 0; i < len; i++)
        {
            arf_srcptr ai, bi, ci, di;
            mp_size_t an, bn, cn, dn;
            slong aexp, bexp, cexp, dexp;

            ai = ((arf_srcptr) x) + 2 * i * xstep;
            bi = ((arf_srcptr) x) + 2 * i * xstep + 1;
            ci = ((arf_srcptr) y) + 2 * i * ystep;
            di = ((arf_srcptr) y) + 2 * i * ystep + 1;

            an = ARF_SIZE(ai);
            bn = ARF_SIZE(bi);
            cn = ARF_SIZE(ci);
            dn = ARF_SIZE(di);

            aexp = ARF_EXP(ai);
            bexp = ARF_EXP(bi);
            cexp = ARF_EXP(ci);
            dexp = ARF_EXP(di);

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
                _arf_complex_mul_gauss(e, f, ai, bi, ci, di);
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
            slong sum_exp;
            mp_ptr sum;
            mp_size_t sn;
            mp_limb_t serr;
            int flipsign;

            if (xoff == yoff)
            {
                sum_exp = re_sum_exp;
                sum = re_sum;
                sn = re_sn;
                if (re_max_exp == WORD_MIN)
                    continue;
            }
            else
            {
                sum_exp = im_sum_exp;
                sum = im_sum;
                sn = im_sn;
                if (im_max_exp == WORD_MIN)
                    continue;
            }

            serr = 0;
            flipsign = (xoff + yoff == 2);

            for (i = 0; i < len; i++)
            {
                xi = ((arf_srcptr) x) + 2 * i * xstep + xoff;
                yi = ((arf_srcptr) y) + 2 * i * ystep + yoff;

                xm = xi;
                ym = yi;

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
                    }
                    else if (xn <= 2 && yn <= 2 && sn <= 3)
                    {
                        mp_limb_t x1, x0, y1, y0;
                        mp_limb_t u3, u2, u1, u0;

                        if (xn == 1 && yn == 1)
                        {
                            x0 = ARF_NOPTR_D(xm)[0];
                            y0 = ARF_NOPTR_D(ym)[0];
                            umul_ppmm(u3, u2, x0, y0);
                            u1 = u0 = 0;
                        }
                        else if (xn == 2 && yn == 2)
                        {
                            x0 = ARF_NOPTR_D(xm)[0];
                            x1 = ARF_NOPTR_D(xm)[1];
                            y0 = ARF_NOPTR_D(ym)[0];
                            y1 = ARF_NOPTR_D(ym)[1];
                            FLINT_MPN_MUL_2X2(u3, u2, u1, u0, x1, x0, y1, y0);
                        }
                        else if (xn == 1)
                        {
                            x0 = ARF_NOPTR_D(xm)[0];
                            y0 = ARF_NOPTR_D(ym)[0];
                            y1 = ARF_NOPTR_D(ym)[1];
                            FLINT_MPN_MUL_2X1(u3, u2, u1, y1, y0, x0);
                            u0 = 0;
                        }
                        else
                        {
                            x0 = ARF_NOPTR_D(xm)[0];
                            x1 = ARF_NOPTR_D(xm)[1];
                            y0 = ARF_NOPTR_D(ym)[0];
                            FLINT_MPN_MUL_2X1(u3, u2, u1, x1, x0, y0);
                            u0 = 0;
                        }

                        if (sn == 2)
                        {
                            if (shift < FLINT_BITS)
                            {
                                u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                                u3 = (u3 >> shift);
                            }
                            else if (shift == FLINT_BITS)
                            {
                                u2 = u3;
                                u3 = 0;
                            }
                            else /* FLINT_BITS < shift < 2 * FLINT_BITS */
                            {
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
                                u1 = (u1 >> shift) | (u2 << (FLINT_BITS - shift));
                                u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                                u3 = (u3 >> shift);
                            }
                            else if (shift == FLINT_BITS)
                            {
                                u1 = u2;
                                u2 = u3;
                                u3 = 0;
                            }
                            else if (shift < 2 * FLINT_BITS)
                            {
                                u1 = (u3 << (2 * FLINT_BITS - shift)) | (u2 >> (shift - FLINT_BITS));
                                u2 = (u3 >> (shift - FLINT_BITS));
                                u3 = 0;
                            }
                            else if (shift == 2 * FLINT_BITS)
                            {
                                u1 = u3;
                                u2 = 0;
                                u3 = 0;
                            }
                            else  /* 2 * FLINT_BITS < shift < 3 * FLINT_BITS */
                            {
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

                        if (use_gauss == NULL || use_gauss[i] == 0)
                            _arb_dot_addmul_generic(sum, &serr, tmp, sn, xptr, xn, yptr, yn, xnegative ^ ynegative ^ flipsign, shift);
                    }
                }
            }
        }
    }

    _arb_dot_output(acf_realref(res), re_sum, re_sn, subtract, re_sum_exp, re_prec, rnd);
    _arb_dot_output(acf_imagref(res), im_sum, im_sn, subtract, im_sum_exp, im_prec, rnd);

    ARF_ADD_TMP_FREE(re_sum, alloc);
    if (use_gauss != NULL)
        flint_free(use_gauss);
}
