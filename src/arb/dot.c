/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mpn_extras.h"

/* We need uint64_t instead of mp_limb_t on 32-bit systems for
   safe summation of 30-bit error bounds. */
#include <stdint.h>

void mpfr_mulhigh_n(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n);
void mpfr_sqrhigh_n(mp_ptr rp, mp_srcptr np, mp_size_t n);

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

/* mag_set_ui_2exp_si but assumes no promotion can occur.
   Do we need this? */
void
mag_set_ui_2exp_small(mag_t z, ulong x, slong e)
{
    _fmpz_demote(MAG_EXPREF(z));

    if (x == 0)
    {
        MAG_EXP(z) = 0;
        MAG_MAN(z) = 0;
    }
    else
    {
        slong bits;
        mp_limb_t overflow;

        bits = flint_clz(x);
        bits = FLINT_BITS - bits;

        if (bits <= MAG_BITS)
        {
            x = x << (MAG_BITS - bits);
        }
        else
        {
            x = (x >> (bits - MAG_BITS)) + 1;
            overflow = x >> MAG_BITS;
            bits += overflow;
            x >>= overflow;
        }

        MAG_EXP(z) = bits + e;
        MAG_MAN(z) = x;
    }
}

/* Sets rad to (Aerr 2^Aexp + Berr 2^Bexp + Cerr 2^Cexp) 2^(-MAG_BITS).
   Assumes that overflow cannot occur. */
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

static void
mulhigh(mp_ptr res, mp_srcptr xptr, mp_size_t xn, mp_srcptr yptr, mp_size_t yn, mp_size_t nn)
{
    mp_ptr tmp, xxx, yyy;
    slong k;
    ARF_MUL_TMP_DECL;

    ARF_MUL_TMP_ALLOC(tmp, 2 * nn);
    xxx = tmp;
    yyy = tmp + nn;

    mpn_copyi(xxx + nn - FLINT_MIN(xn, nn), xptr + xn - FLINT_MIN(xn, nn), FLINT_MIN(xn, nn));
    mpn_copyi(yyy + nn - FLINT_MIN(yn, nn), yptr + yn - FLINT_MIN(yn, nn), FLINT_MIN(yn, nn));

    for (k = 0; k < nn - xn; k++)
        xxx[k] = 0;
    for (k = 0; k < nn - yn; k++)
        yyy[k] = 0;

    if (xptr == yptr && xn == yn)
        mpfr_sqrhigh_n(res, xxx, nn);
    else
        mpfr_mulhigh_n(res, xxx, yyy, nn);

    ARF_MUL_TMP_FREE(tmp, 2 * nn);
}

void
_arb_dot_addmul_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn, mp_srcptr yptr, mp_size_t yn,
    int negative, flint_bitcnt_t shift)
{
    slong shift_bits, shift_limbs, term_prec;
    mp_limb_t cy;
    mp_ptr sstart, tstart;
    mp_size_t tn, nn;

    shift_bits = shift % FLINT_BITS;
    shift_limbs = shift / FLINT_BITS;

    /*
            shift          term_prec      (discarded bits)
    |--------------------|--------------|
                         | t[tn-1] |   ...  |  t[0]  |         Term
    [ s[n-1] | s[n-2] |   ...  |  s[0]  ]                      Sum

    */

    /* Upper bound for the term bit length. The actual term bit length
       could be smaller if the product does not extend all
       the way down to s. */
    term_prec = sn * FLINT_BITS - shift;

    /* The mulhigh error relative to the top of the sum will
       be bounded by nn * 2^(-shift) * 2^(-FLINT_BITS*nn).
       One extra limb ensures that this is <1 sum-ulp. */
    term_prec += FLINT_BITS;
    nn = (term_prec + FLINT_BITS - 1) / FLINT_BITS;

    /* Sanity check; must conform to the pre-allocated memory! */
    if (nn > sn + 2)
    {
        flint_throw(FLINT_ERROR, "nn > sn + 2\n");
    }

    /* Use mulhigh? */
    if (term_prec >= MUL_MPFR_MIN_LIMBS * FLINT_BITS &&
        term_prec <= MUL_MPFR_MAX_LIMBS * FLINT_BITS &&
        xn * FLINT_BITS > 0.9 * term_prec &&
        yn * FLINT_BITS > 0.9 * term_prec)
    {
        mulhigh(tmp, xptr, xn, yptr, yn, nn);
        tstart = tmp + nn;
        tn = nn;
        serr[0]++;
    }
    else
    {
        if (xn > nn || yn > nn)
        {
            if (xn > nn)
            {
                xptr += (xn - nn);
                xn = nn;
            }

            if (yn > nn)
            {
                yptr += (yn - nn);
                yn = nn;
            }

            serr[0]++;
        }

        tn = xn + yn;
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp + 1, xptr, xn, yptr, yn);
        tstart = tmp + 1;
    }

    if (shift_bits != 0)
    {
        tstart[-1] = mpn_rshift(tstart, tstart, tn, shift_bits);
        tstart = tstart - 1;
        tn++;
    }

    while (tstart[0] == 0)
    {
        tstart++;
        tn--;
    }

    if (shift_limbs + tn <= sn)
    {
        /* No truncation of the term. */
        sstart = sum + sn - shift_limbs - tn;
        nn = tn;
    }
    else
    {
        /* The term is truncated. */
        sstart = sum;
        tstart = tstart - (sn - shift_limbs - tn);
        nn = sn - shift_limbs;
        serr[0]++;
    }

    /* Likely case. Note: assumes 1 extra (dummy) limb in sum. */
    if (shift_limbs <= 1)
    {
        if (negative)
            sstart[nn] -= mpn_sub_n(sstart, sstart, tstart, nn);
        else
            sstart[nn] += mpn_add_n(sstart, sstart, tstart, nn);
    }
    else
    {

        if (negative)
        {
            cy = mpn_sub_n(sstart, sstart, tstart, nn);
            mpn_sub_1(sstart + nn, sstart + nn, shift_limbs, cy);
        }
        else
        {
            cy = mpn_add_n(sstart, sstart, tstart, nn);
            mpn_add_1(sstart + nn, sstart + nn, shift_limbs, cy);
        }
    }
}

void
_arb_dot_add_generic(mp_ptr sum, mp_ptr serr, mp_ptr tmp, mp_size_t sn,
    mp_srcptr xptr, mp_size_t xn,
    int negative, flint_bitcnt_t shift)
{
    slong shift_bits, shift_limbs, term_prec;
    mp_limb_t cy, err;
    mp_ptr sstart, tstart;
    mp_size_t tn, nn;

    shift_bits = shift % FLINT_BITS;
    shift_limbs = shift / FLINT_BITS;

    term_prec = sn * FLINT_BITS - shift;
    term_prec += FLINT_BITS;
    nn = (term_prec + FLINT_BITS - 1) / FLINT_BITS;

    err = 0;

    if (xn > nn)
    {
        xptr += (xn - nn);
        xn = nn;
        err = 1;
    }

    tn = xn;

    if (shift_bits == 0)
    {
        mpn_copyi(tmp, xptr, tn);
        tstart = tmp;
    }
    else
    {
        tmp[0] = mpn_rshift(tmp + 1, xptr, tn, shift_bits);
        tstart = tmp;
        tn++;
    }

    while (tstart[0] == 0)
    {
        tstart++;
        tn--;
    }

    if (shift_limbs + tn <= sn)
    {
        /* No truncation of the term. */
        sstart = sum + sn - shift_limbs - tn;
        nn = tn;
    }
    else
    {
        /* The term is truncated. */
        sstart = sum;
        tstart = tstart - (sn - shift_limbs - tn);
        nn = sn - shift_limbs;
        err = 1;
    }

    serr[0] += err;

    /* Likely case. Note: assumes 1 extra (dummy) limb in sum. */
    if (shift_limbs <= 1)
    {
        if (negative)
            sstart[nn] -= mpn_sub_n(sstart, sstart, tstart, nn);
        else
            sstart[nn] += mpn_add_n(sstart, sstart, tstart, nn);
    }
    else
    {

        if (negative)
        {
            cy = mpn_sub_n(sstart, sstart, tstart, nn);
            mpn_sub_1(sstart + nn, sstart + nn, shift_limbs, cy);
        }
        else
        {
            cy = mpn_add_n(sstart, sstart, tstart, nn);
            mpn_add_1(sstart + nn, sstart + nn, shift_limbs, cy);
        }
    }
}

void
arb_dot(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
{
    slong i, j, nonzero, padding, extend;
    slong xexp, yexp, exp, max_exp, min_exp, sum_exp;
    slong xrexp, yrexp, srad_exp, max_rad_exp;
    int xnegative, ynegative, inexact;
    mp_size_t xn, yn, sn, alloc;
    flint_bitcnt_t shift;
    arb_srcptr xi, yi;
    arf_srcptr xm, ym;
    mag_srcptr xr, yr;
    mp_limb_t xtop, ytop;
    mp_limb_t xrad, yrad;
    mp_limb_t serr;   /* Sum over arithmetic errors */
    uint64_t srad;    /* Sum over propagated errors */
    mp_ptr tmp, sum;  /* Workspace */
    ARF_ADD_TMP_DECL;

    /* todo: fast fma and fmma (len=2) code */
    if (len <= 1)
    {
        if (initial == NULL)
        {
            if (len <= 0)
                arb_zero(res);
            else
            {
                arb_mul(res, x, y, prec);
                if (subtract)
                    arb_neg(res, res);
            }
            return;
        }
        else if (len <= 0)
        {
            arb_set_round(res, initial, prec);
            return;
        }
    }

    /* Number of nonzero midpoint terms in sum. */
    nonzero = 0;

    /* Terms are bounded by 2^max_exp (with WORD_MIN = -infty) */
    max_exp = WORD_MIN;

    /* Propagated error terms are bounded by 2^max_rad_exp */
    max_rad_exp = WORD_MIN;

    /* Used to reduce the precision. */
    min_exp = WORD_MAX;

    /* Account for the initial term. */
    if (initial != NULL)
    {
        if (!ARB_IS_LAGOM(initial))
        {
            arb_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec);
            return;
        }

        xm = arb_midref(initial);
        xr = arb_radref(initial);

        if (!arf_is_special(xm))
        {
            max_exp = ARF_EXP(xm);
            nonzero++;

            if (prec > 2 * FLINT_BITS)
                min_exp = ARF_EXP(xm) - ARF_SIZE(xm) * FLINT_BITS;
        }

        if (!mag_is_special(xr))
            max_rad_exp = MAG_EXP(xr);
    }

    /* Determine maximum exponents for the main sum and the radius sum. */
    for (i = 0; i < len; i++)
    {
        xi = x + i * xstep;
        yi = y + i * ystep;

        /* Fallback for huge exponents or non-finite values. */
        if (!ARB_IS_LAGOM(xi) || !ARB_IS_LAGOM(yi))
        {
            arb_dot_simple(res, initial, subtract, x, xstep, y, ystep, len, prec);
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

    /* The midpoint sum is zero. */
    if (max_exp == WORD_MIN)
    {
        /* The sum is exactly zero. */
        if (max_rad_exp == WORD_MIN)
        {
            arb_zero(res);
            return;
        }

        prec = 2;
    }
    else
    {
        /* Reduce precision based on errors. */
        if (max_rad_exp != WORD_MIN)
            prec = FLINT_MIN(prec, max_exp - max_rad_exp + MAG_BITS);

        /* Reduce precision based on actual sizes. */
        if (min_exp != WORD_MAX)
            prec = FLINT_MIN(prec, max_exp - min_exp + MAG_BITS);

        prec = FLINT_MAX(prec, 2);
    }

    /* Extend sum so that we can use two's complement addition. */
    extend = FLINT_BIT_COUNT(nonzero) + 1;

    /* Extra bits to improve accuracy (optional). */
    padding = 4 + FLINT_BIT_COUNT(len);

    /* Number of limbs. */
    sn = (prec + extend + padding + FLINT_BITS - 1) / FLINT_BITS;

    /* Avoid having to make a special case for sn = 1. */
    sn = FLINT_MAX(sn, 2);

    /* Exponent for the main sum. */
    sum_exp = max_exp + extend;

    /* We need sn + 1 limb for the sum (sn limbs + 1 dummy limb
       for carry or borrow that avoids an extra branch). We need
       2 * (sn + 2) limbs to store the product of two numbers
       with up to (sn + 2) limbs, plus 1 extra limb for shifting
       the product. */
    alloc = (sn + 1) + 2 * (sn + 2) + 1;
    ARF_ADD_TMP_ALLOC(sum, alloc)
    tmp = sum + (sn + 1);

    /* Sum of propagated errors. */
    srad_exp = max_rad_exp;
    srad = 0;

    /* Set sum to 0 */
    serr = 0;
    for (j = 0; j < sn + 1; j++)
        sum[j] = 0;

    if (initial != NULL)
    {
        xm = arb_midref(initial);
        xr = arb_radref(initial);

        if (!arf_is_special(xm))
        {
            mp_srcptr xptr;

            xexp = ARF_EXP(xm);
            xn = ARF_SIZE(xm);
            xnegative = ARF_SGNBIT(xm);

            shift = sum_exp - xexp;

            if (shift >= sn * FLINT_BITS)
            {
                serr++;
            }
            else
            {
                xptr = (xn <= ARF_NOPTR_LIMBS) ? ARF_NOPTR_D(xm) : ARF_PTR_D(xm);
                _arb_dot_add_generic(sum, &serr, tmp, sn, xptr, xn, xnegative ^ subtract, shift);
            }
        }

        if (!mag_is_special(xr))
        {
            xrad = MAG_MAN(xr);
            xrexp = MAG_EXP(xr);

            shift = srad_exp - xrexp;
            if (shift < 64)
                srad += (xrad >> shift) + 1;
            else
                srad++;
        }
    }

    for (i = 0; i < len; i++)
    {
        xi = x + i * xstep;
        yi = y + i * ystep;
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
#if 0
            else if (xn == 1 && yn == 1 && sn == 2 && shift < FLINT_BITS)  /* Fastest path. */
            {
                mp_limb_t hi, lo, out;

                xtop = ARF_NOPTR_D(xm)[0];
                ytop = ARF_NOPTR_D(ym)[0];

                umul_ppmm(hi, lo, xtop, ytop);

                out = lo << (FLINT_BITS - shift);
                lo = (lo >> shift) | (hi << (FLINT_BITS - shift));
                hi = (hi >> shift);
                serr += (out != 0);

                if (xnegative ^ ynegative)
                    sub_ddmmss(sum[1], sum[0], sum[1], sum[0], hi, lo);
                else
                    add_ssaaaa(sum[1], sum[0], sum[1], sum[0], hi, lo);
            }
            else if (xn == 2 && yn == 2 && shift < FLINT_BITS && sn <= 3)
            {
                mp_limb_t x1, x0, y1, y0;
                mp_limb_t u3, u2, u1, u0;

                x0 = ARF_NOPTR_D(xm)[0];
                x1 = ARF_NOPTR_D(xm)[1];
                y0 = ARF_NOPTR_D(ym)[0];
                y1 = ARF_NOPTR_D(ym)[1];

                xtop = x1;
                ytop = y1;

                FLINT_MPN_MUL_2X2(u3, u2, u1, u0, x1, x0, y1, y0);

                u0 = (u0 != 0) || ((u1 << (FLINT_BITS - shift)) != 0);
                u1 = (u1 >> shift) | (u2 << (FLINT_BITS - shift));
                u2 = (u2 >> shift) | (u3 << (FLINT_BITS - shift));
                u3 = (u3 >> shift);

                if (sn == 2)
                {
                    serr += (u0 || (u1 != 0));
                    if (xnegative ^ ynegative)
                        sub_ddmmss(sum[1], sum[0], sum[1], sum[0], u3, u2);
                    else
                        add_ssaaaa(sum[1], sum[0], sum[1], sum[0], u3, u2);
                }
                else
                {
                    serr += u0;
                    if (xnegative ^ ynegative)
                        sub_dddmmmsss(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                    else
                        add_sssaaaaaa(sum[2], sum[1], sum[0], sum[2], sum[1], sum[0], u3, u2, u1);
                }
            }
#endif
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

                    if (xnegative ^ ynegative)
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

                    if (xnegative ^ ynegative)
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

                _arb_dot_addmul_generic(sum, &serr, tmp, sn, xptr, xn, yptr, yn, xnegative ^ ynegative, shift);
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

    xnegative = 0;
    if (sum[sn - 1] >= LIMB_TOP)
    {
        mpn_neg(sum, sum, sn);
        xnegative = 1;
    }

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
            inexact = _arf_set_round_mpn(arb_midref(res), &exp, sum, sn2, xnegative ^ subtract, prec, ARF_RND_DOWN);
            _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp + sum_exp2);
        }
    }
    else
    {

        if (sn == 2)
            inexact = _arf_set_round_uiui(arb_midref(res), &exp, sum[1], sum[0], xnegative ^ subtract, prec, ARF_RND_DOWN);
        else
            inexact = _arf_set_round_mpn(arb_midref(res), &exp, sum, sn, xnegative ^ subtract, prec, ARF_RND_DOWN);

        _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp + sum_exp);
    }

    /*
    Add the three sources of error.

    Final rounding error   = inexact * 2^(exp + sum_exp - prec)
                             0 <= inexact <= 1

    Arithmetic error       = serr    * 2^(sum_exp - sn * FLINT_BITS)
                             0 <= serr <= 3 * n

    Propagated error       = srad    * 2^(srad_exp - MAG_BITS)
                             0 <= srad <= 3 * n * MAG_BITS

    We shift the first two by MAG_BITS so that the magnitudes
    become similar and a reasonably accurate addition can be done
    by comparing exponents without normalizing first.
    */

    add_errors(arb_radref(res),
        inexact << MAG_BITS,
        exp + sum_exp - prec,
        ((uint64_t) serr) << MAG_BITS,
        sum_exp - sn * FLINT_BITS,
        srad,
        srad_exp);

    ARF_ADD_TMP_FREE(sum, alloc);
}
