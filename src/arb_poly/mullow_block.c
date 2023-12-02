/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "arb_poly.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
#else
# include <math.h>
#endif

void
_arb_poly_get_scale(fmpz_t scale, arb_srcptr x, slong xlen,
                                  arb_srcptr y, slong ylen)
{
    slong xa, xb, ya, yb, den;

    fmpz_zero(scale);

    /* ignore zeros (and infs/nans!); find the first and last
       finite nonzero entries to determine the scale */
    xa = 0;
    xb = xlen - 1;
    while (xa < xlen && arf_is_special(arb_midref(x + xa))) xa++;
    while (xb > xa && arf_is_special(arb_midref(x + xb))) xb--;

    ya = 0;
    yb = ylen - 1;
    while (ya < ylen && arf_is_special(arb_midref(y + ya))) ya++;
    while (yb > ya && arf_is_special(arb_midref(y + yb))) yb--;

    /* compute average of exponent differences, weighted by the lengths */
    if (xa <= xb && ya <= yb && (xa < xb || ya < yb))
    {
        fmpz_add(scale, scale, ARF_EXPREF(arb_midref(x + xb)));
        fmpz_sub(scale, scale, ARF_EXPREF(arb_midref(x + xa)));
        fmpz_add(scale, scale, ARF_EXPREF(arb_midref(y + yb)));
        fmpz_sub(scale, scale, ARF_EXPREF(arb_midref(y + ya)));

        den = (xb - xa) + (yb - ya);

        /* scale = floor(scale / den + 1/2) = floor((2 scale + den) / (2 den)) */
        fmpz_mul_2exp(scale, scale, 1);
        fmpz_add_ui(scale, scale, den);
        fmpz_fdiv_q_ui(scale, scale, 2 * den);
    }
}

/* Break vector into same-exponent blocks where the largest block
   has a height of at most ALPHA*prec + BETA bits. These are just
   tuning parameters. Note that ALPHA * MAG_BITS + BETA
   should be smaller than DOUBLE_BLOCK_MAX_HEIGHT if we want to use
   doubles for error bounding. */
#define ALPHA 3.0
#define BETA 512


/* Maximum length of block for which we use double multiplication
   (for longer blocks, we use fmpz_poly multiplication). This is essentially
   just a tuning parameter, but note that it must be considered when
   compensating for rounding error below. */
#define DOUBLE_BLOCK_MAX_LENGTH 1000

/* Computing a dot product of length DOUBLE_BLOCK_MAX_LENGTH involving
   only nonnegative numbers, and then multiplying by this factor, must give
   an upper bound for the exact dot product (we can assume that no
   overflow or underflow occurs). The following is certainly
   sufficient, but it would be nice to include a formal proof here. */
#define DOUBLE_ROUNDING_FACTOR (1.0 + 1e-9)

/* Maximum height for which we use double multiplication. Since the dynamic
   exponent range of doubles is about +/- 1024, this must be less than about
   1024 (to allow the product of two numbers). This must also
   account for adding MAG_BITS bits. */
#define DOUBLE_BLOCK_MAX_HEIGHT 800

/* We divide coefficients by 2^DOUBLE_BLOCK_SHIFT when converting them to
   doubles, in order to use the whole exponent range. Note that this means
   numbers of size (2^(-DOUBLE_BLOCK_SHIFT))^2 must not underflow. */
#define DOUBLE_BLOCK_SHIFT (DOUBLE_BLOCK_MAX_HEIGHT / 2)


static void
_mag_vec_get_fmpz_2exp_blocks(fmpz * coeffs,
    double * dblcoeffs, fmpz * exps, slong * blocks, const fmpz_t scale,
    arb_srcptr x, mag_srcptr xm, slong len)
{
    fmpz_t top, bot, t, b, v, block_top, block_bot;
    slong i, j, s, block, bits, maxheight;
    int in_zero;
    mag_srcptr cur;

    fmpz_init(top);
    fmpz_init(bot);
    fmpz_init(t);
    fmpz_init(b);
    fmpz_init(v);
    fmpz_init(block_top);
    fmpz_init(block_bot);

    blocks[0] = 0;
    block = 0;
    in_zero = 1;

    maxheight = ALPHA * MAG_BITS + BETA;
    if (maxheight > DOUBLE_BLOCK_MAX_HEIGHT)
        flint_throw(FLINT_ERROR, "(%s): maxheight > DOUBLE_BLOCK_MAX_HEIGHT\n", __func__);

    for (i = 0; i < len; i++)
    {
        cur = (x == NULL) ? (xm + i) : arb_radref(x + i);

        /* Skip (must be zero, since we assume there are no Infs/NaNs). */
        if (mag_is_special(cur))
            continue;

        /* Bottom and top exponent of current number */
        bits = MAG_BITS;
        fmpz_set(top, MAG_EXPREF(cur));
        fmpz_submul_ui(top, scale, i);
        fmpz_sub_ui(bot, top, bits);

        /* Extend current block. */
        if (in_zero)
        {
            fmpz_swap(block_top, top);
            fmpz_swap(block_bot, bot);
        }
        else
        {
            fmpz_max(t, top, block_top);
            fmpz_min(b, bot, block_bot);
            fmpz_sub(v, t, b);

            /* extend current block */
            if (fmpz_cmp_ui(v, maxheight) < 0)
            {
                fmpz_swap(block_top, t);
                fmpz_swap(block_bot, b);
            }
            else  /* start new block */
            {
                /* write exponent for previous block */
                fmpz_set(exps + block, block_bot);

                block++;
                blocks[block] = i;

                fmpz_swap(block_top, top);
                fmpz_swap(block_bot, bot);
            }
        }

        in_zero = 0;
    }

    /* write exponent for last block */
    fmpz_set(exps + block, block_bot);

    /* end marker */
    blocks[block + 1] = len;

    /* write the block data */
    for (i = 0; blocks[i] != len; i++)
    {
        for (j = blocks[i]; j < blocks[i + 1]; j++)
        {
            cur = (x == NULL) ? (xm + j) : arb_radref(x + j);

            if (mag_is_special(cur))
            {
                fmpz_zero(coeffs + j);
                dblcoeffs[j] = 0.0;
            }
            else
            {
                mp_limb_t man;
                double c;

                man = MAG_MAN(cur);

                /* TODO: only write and use doubles when block is short? */

                /* Divide by 2^(scale * j) */
                fmpz_mul_ui(t, scale, j);
                fmpz_sub(t, MAG_EXPREF(cur), t);

                fmpz_sub_ui(t, t, MAG_BITS); /* bottom exponent */
                s = _fmpz_sub_small(t, exps + i);

                if (s < 0)
                    flint_throw(FLINT_ERROR, "(%s): s < 0\n", __func__);

                fmpz_set_ui(coeffs + j, man);
                fmpz_mul_2exp(coeffs + j, coeffs + j, s);
                c = man;
                c = ldexp(c, s - DOUBLE_BLOCK_SHIFT);
                if (c < 1e-150 || c > 1e150) /* Bug catcher */
                    flint_throw(FLINT_ERROR, "(%s): c large or big\n", __func__);
                dblcoeffs[j] = c;
            }
        }
    }

    fmpz_clear(top);
    fmpz_clear(bot);
    fmpz_clear(t);
    fmpz_clear(b);
    fmpz_clear(v);
    fmpz_clear(block_top);
    fmpz_clear(block_bot);
}

static void
_arb_vec_get_fmpz_2exp_blocks(fmpz * coeffs, fmpz * exps,
    slong * blocks, const fmpz_t scale, arb_srcptr x, slong len, slong prec)
{
    fmpz_t top, bot, t, b, v, block_top, block_bot;
    slong i, j, s, block, bits, maxheight;
    int in_zero;

    fmpz_init(top);
    fmpz_init(bot);
    fmpz_init(t);
    fmpz_init(b);
    fmpz_init(v);
    fmpz_init(block_top);
    fmpz_init(block_bot);

    blocks[0] = 0;
    block = 0;
    in_zero = 1;

    if (prec == ARF_PREC_EXACT)
        maxheight = ARF_PREC_EXACT;
    else
        maxheight = ALPHA * prec + BETA;

    for (i = 0; i < len; i++)
    {
        bits = arf_bits(arb_midref(x + i));

        /* Skip (must be zero, since we assume there are no Infs/NaNs). */
        if (bits == 0)
            continue;

        /* Bottom and top exponent of current number */
        fmpz_set(top, ARF_EXPREF(arb_midref(x + i)));
        fmpz_submul_ui(top, scale, i);
        fmpz_sub_ui(bot, top, bits);

        /* Extend current block. */
        if (in_zero)
        {
            fmpz_swap(block_top, top);
            fmpz_swap(block_bot, bot);
        }
        else
        {
            fmpz_max(t, top, block_top);
            fmpz_min(b, bot, block_bot);
            fmpz_sub(v, t, b);

            /* extend current block */
            if (fmpz_cmp_ui(v, maxheight) < 0)
            {
                fmpz_swap(block_top, t);
                fmpz_swap(block_bot, b);
            }
            else  /* start new block */
            {
                /* write exponent for previous block */
                fmpz_set(exps + block, block_bot);

                block++;
                blocks[block] = i;

                fmpz_swap(block_top, top);
                fmpz_swap(block_bot, bot);
            }
        }

        in_zero = 0;
    }

    /* write exponent for last block */
    fmpz_set(exps + block, block_bot);

    /* end marker */
    blocks[block + 1] = len;

    /* write the block data */
    for (i = 0; blocks[i] != len; i++)
    {
        for (j = blocks[i]; j < blocks[i + 1]; j++)
        {
            if (arf_is_special(arb_midref(x + j)))
            {
                fmpz_zero(coeffs + j);
            }
            else
            {
                /* TODO: make this a single operation */
                arf_get_fmpz_2exp(coeffs + j, bot, arb_midref(x + j));

                fmpz_mul_ui(t, scale, j);
                fmpz_sub(t, bot, t);
                s = _fmpz_sub_small(t, exps + i);
                if (s < 0)
                    flint_throw(FLINT_ERROR, "(%s): s < 0\n", __func__);
                fmpz_mul_2exp(coeffs + j, coeffs + j, s);
            }
        }
    }

    fmpz_clear(top);
    fmpz_clear(bot);
    fmpz_clear(t);
    fmpz_clear(b);
    fmpz_clear(v);
    fmpz_clear(block_top);
    fmpz_clear(block_bot);
}

static void
_arb_poly_addmullow_rad(arb_ptr z, fmpz * zz,
    const fmpz * xz, const double * xdbl, const fmpz * xexps,
    const slong * xblocks, slong xlen,
    const fmpz * yz, const double * ydbl, const fmpz * yexps,
    const slong * yblocks, slong ylen, slong n)
{
    slong i, j, k, ii, xp, yp, xl, yl, bn;
    fmpz_t zexp;
    mag_t t;

    fmpz_init(zexp);
    mag_init(t);

    for (i = 0; (xp = xblocks[i]) != xlen; i++)
    {
        for (j = 0; (yp = yblocks[j]) != ylen; j++)
        {
            if (xp + yp >= n)
                continue;

            xl = xblocks[i + 1] - xp;
            yl = yblocks[j + 1] - yp;
            bn = FLINT_MIN(xl + yl - 1, n - xp - yp);
            xl = FLINT_MIN(xl, bn);
            yl = FLINT_MIN(yl, bn);

            fmpz_add_inline(zexp, xexps + i, yexps + j);

            if (xl > 1 && yl > 1 &&
                (xl < DOUBLE_BLOCK_MAX_LENGTH || yl < DOUBLE_BLOCK_MAX_LENGTH))
            {
                fmpz_add_ui(zexp, zexp, 2 * DOUBLE_BLOCK_SHIFT);

                for (k = 0; k < bn; k++)
                {
                    /* Classical multiplication (may round down!) */
                    double ss = 0.0;

                    for (ii = FLINT_MAX(0, k - yl + 1);
                        ii <= FLINT_MIN(xl - 1, k); ii++)
                    {
                        ss += xdbl[xp + ii] * ydbl[yp + k - ii];
                    }

                    /* Compensate for rounding error */
                    ss *= DOUBLE_ROUNDING_FACTOR;

                    mag_set_d_2exp_fmpz(t, ss, zexp);
                    mag_add(arb_radref(z + xp + yp + k),
                            arb_radref(z + xp + yp + k), t);
                }
            }
            else
            {
                if (xl >= yl)
                    _fmpz_poly_mullow(zz, xz + xp, xl, yz + yp, yl, bn);
                else
                    _fmpz_poly_mullow(zz, yz + yp, yl, xz + xp, xl, bn);

                for (k = 0; k < bn; k++)
                {
                    mag_set_fmpz_2exp_fmpz(t, zz + k, zexp);
                    mag_add(arb_radref(z + xp + yp + k),
                            arb_radref(z + xp + yp + k), t);
                }
            }
        }
    }

    fmpz_clear(zexp);
    mag_clear(t);
}

static void
_arb_poly_addmullow_block(arb_ptr z, fmpz * zz,
    const fmpz * xz, const fmpz * xexps, const slong * xblocks, slong xlen,
    const fmpz * yz, const fmpz * yexps, const slong * yblocks, slong ylen,
    slong n, slong prec, int squaring)
{
    slong i, j, k, xp, yp, xl, yl, bn;
    fmpz_t zexp;

    fmpz_init(zexp);

    if (squaring)
    {
        for (i = 0; (xp = xblocks[i]) != xlen; i++)
        {
            if (2 * xp >= n)
                continue;

            xl = xblocks[i + 1] - xp;
            bn = FLINT_MIN(2 * xl - 1, n - 2 * xp);
            xl = FLINT_MIN(xl, bn);

            _fmpz_poly_sqrlow(zz, xz + xp, xl, bn);
            _fmpz_add2_fast(zexp, xexps + i, xexps + i, 0);

            for (k = 0; k < bn; k++)
                arb_add_fmpz_2exp(z + 2 * xp + k, z + 2 * xp + k, zz + k, zexp, prec);
        }
    }

    for (i = 0; (xp = xblocks[i]) != xlen; i++)
    {
        for (j = squaring ? i + 1 : 0; (yp = yblocks[j]) != ylen; j++)
        {
            if (xp + yp >= n)
                continue;

            xl = xblocks[i + 1] - xp;
            yl = yblocks[j + 1] - yp;
            bn = FLINT_MIN(xl + yl - 1, n - xp - yp);
            xl = FLINT_MIN(xl, bn);
            yl = FLINT_MIN(yl, bn);

            if (xl >= yl)
                _fmpz_poly_mullow(zz, xz + xp, xl, yz + yp, yl, bn);
            else
                _fmpz_poly_mullow(zz, yz + yp, yl, xz + xp, xl, bn);

           _fmpz_add2_fast(zexp, xexps + i, yexps + j, squaring);

            for (k = 0; k < bn; k++)
                arb_add_fmpz_2exp(z + xp + yp + k, z + xp + yp + k, zz + k, zexp, prec);
        }
    }

    fmpz_clear(zexp);
}

void
_arb_poly_mullow_block(arb_ptr z, arb_srcptr x, slong xlen,
                                arb_srcptr y, slong ylen, slong n, slong prec)
{
    slong xmlen, xrlen, ymlen, yrlen, i;
    fmpz *xz, *yz, *zz;
    fmpz *xe, *ye;
    slong *xblocks, *yblocks;
    int squaring;
    fmpz_t scale, t;

    xlen = FLINT_MIN(xlen, n);
    ylen = FLINT_MIN(ylen, n);

    squaring = (x == y) && (xlen == ylen);

    /* Strip trailing zeros */
    xmlen = xrlen = xlen;
    while (xmlen > 0 && arf_is_zero(arb_midref(x + xmlen - 1))) xmlen--;
    while (xrlen > 0 && mag_is_zero(arb_radref(x + xrlen - 1))) xrlen--;

    if (squaring)
    {
        ymlen = xmlen;
        yrlen = xrlen;
    }
    else
    {
        ymlen = yrlen = ylen;
        while (ymlen > 0 && arf_is_zero(arb_midref(y + ymlen - 1))) ymlen--;
        while (yrlen > 0 && mag_is_zero(arb_radref(y + yrlen - 1))) yrlen--;
    }

    /* We don't know how to deal with infinities or NaNs */
    if (!_arb_vec_is_finite(x, xlen) ||
        (!squaring && !_arb_vec_is_finite(y, ylen)))
    {
        _arb_poly_mullow_classical(z, x, xlen, y, ylen, n, prec);
        return;
    }

    xlen = FLINT_MAX(xmlen, xrlen);
    ylen = FLINT_MAX(ymlen, yrlen);

    /* Start with the zero polynomial */
    _arb_vec_zero(z, n);

    /* Nothing to do */
    if (xlen == 0 || ylen == 0)
        return;

    n = FLINT_MIN(n, xlen + ylen - 1);

    fmpz_init(scale);
    fmpz_init(t);
    xz = _fmpz_vec_init(xlen);
    yz = _fmpz_vec_init(ylen);
    zz = _fmpz_vec_init(n);
    xe = _fmpz_vec_init(xlen);
    ye = _fmpz_vec_init(ylen);
    xblocks = flint_malloc(sizeof(slong) * (xlen + 1));
    yblocks = flint_malloc(sizeof(slong) * (ylen + 1));

    _arb_poly_get_scale(scale, x, xlen, y, ylen);

    /* Error propagation */
    /* (xm + xr)*(ym + yr) = (xm*ym) + (xr*ym + xm*yr + xr*yr)
                           = (xm*ym) + (xm*yr + xr*(ym + yr))  */
    if (xrlen != 0 || yrlen != 0)
    {
        mag_ptr tmp;
        double *xdbl, *ydbl;

        tmp = _mag_vec_init(FLINT_MAX(xlen, ylen));
        xdbl = flint_malloc(sizeof(double) * xlen);
        ydbl = flint_malloc(sizeof(double) * ylen);

        /* (xm + xr)^2 = (xm*ym) + (xr^2 + 2 xm xr)
                       = (xm*ym) + xr*(2 xm + xr)    */
        if (squaring)
        {
            _mag_vec_get_fmpz_2exp_blocks(xz, xdbl, xe, xblocks, scale, x, NULL, xrlen);

            for (i = 0; i < xlen; i++)
            {
                arf_get_mag(tmp + i, arb_midref(x + i));
                mag_mul_2exp_si(tmp + i, tmp + i, 1);
                mag_add(tmp + i, tmp + i, arb_radref(x + i));
            }

            _mag_vec_get_fmpz_2exp_blocks(yz, ydbl, ye, yblocks, scale, NULL, tmp, xlen);
            _arb_poly_addmullow_rad(z, zz, xz, xdbl, xe, xblocks, xrlen, yz, ydbl, ye, yblocks, xlen, n);
        }
        else if (yrlen == 0)
        {
            /* xr * |ym| */
            _mag_vec_get_fmpz_2exp_blocks(xz, xdbl, xe, xblocks, scale, x, NULL, xrlen);

            for (i = 0; i < ymlen; i++)
                arf_get_mag(tmp + i, arb_midref(y + i));

            _mag_vec_get_fmpz_2exp_blocks(yz, ydbl, ye, yblocks, scale, NULL, tmp, ymlen);
            _arb_poly_addmullow_rad(z, zz, xz, xdbl, xe, xblocks, xrlen, yz, ydbl, ye, yblocks, ymlen, n);
        }
        else
        {
            /* |xm| * yr */
            for (i = 0; i < xmlen; i++)
                arf_get_mag(tmp + i, arb_midref(x + i));

            _mag_vec_get_fmpz_2exp_blocks(xz, xdbl, xe, xblocks, scale, NULL, tmp, xmlen);
            _mag_vec_get_fmpz_2exp_blocks(yz, ydbl, ye, yblocks, scale, y, NULL, yrlen);
            _arb_poly_addmullow_rad(z, zz, xz, xdbl, xe, xblocks, xmlen, yz, ydbl, ye, yblocks, yrlen, n);

            /* xr*(|ym| + yr) */
            if (xrlen != 0)
            {
                _mag_vec_get_fmpz_2exp_blocks(xz, xdbl, xe, xblocks, scale, x, NULL, xrlen);

                for (i = 0; i < ylen; i++)
                    arb_get_mag(tmp + i, y + i);

                _mag_vec_get_fmpz_2exp_blocks(yz, ydbl, ye, yblocks, scale, NULL, tmp, ylen);
                _arb_poly_addmullow_rad(z, zz, xz, xdbl, xe, xblocks, xrlen, yz, ydbl, ye, yblocks, ylen, n);
            }
        }

        _mag_vec_clear(tmp, FLINT_MAX(xlen, ylen));
        flint_free(xdbl);
        flint_free(ydbl);
    }

    /* multiply midpoints */
    if (xmlen != 0 && ymlen != 0)
    {
        _arb_vec_get_fmpz_2exp_blocks(xz, xe, xblocks, scale, x, xmlen, prec);

        if (squaring)
        {
            _arb_poly_addmullow_block(z, zz, xz, xe, xblocks, xmlen, xz, xe, xblocks, xmlen, n, prec, 1);
        }
        else
        {
            _arb_vec_get_fmpz_2exp_blocks(yz, ye, yblocks, scale, y, ymlen, prec);
            _arb_poly_addmullow_block(z, zz, xz, xe, xblocks, xmlen, yz, ye, yblocks, ymlen, n, prec, 0);
        }
    }

    /* Unscale. */
    if (!fmpz_is_zero(scale))
    {
        fmpz_zero(t);
        for (i = 0; i < n; i++)
        {
            arb_mul_2exp_fmpz(z + i, z + i, t);
            fmpz_add(t, t, scale);
        }
    }

    _fmpz_vec_clear(xz, xlen);
    _fmpz_vec_clear(yz, ylen);
    _fmpz_vec_clear(zz, n);
    _fmpz_vec_clear(xe, xlen);
    _fmpz_vec_clear(ye, ylen);
    flint_free(xblocks);
    flint_free(yblocks);
    fmpz_clear(scale);
    fmpz_clear(t);
}

void
arb_poly_mullow_block(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong n, slong prec)
{
    slong xlen, ylen, zlen;

    xlen = poly1->length;
    ylen = poly2->length;

    if (xlen == 0 || ylen == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    xlen = FLINT_MIN(xlen, n);
    ylen = FLINT_MIN(ylen, n);
    zlen = FLINT_MIN(xlen + ylen - 1, n);

    if (res == poly1 || res == poly2)
    {
        arb_poly_t tmp;
        arb_poly_init2(tmp, zlen);
        _arb_poly_mullow_block(tmp->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
        arb_poly_swap(res, tmp);
        arb_poly_clear(tmp);
    }
    else
    {
        arb_poly_fit_length(res, zlen);
        _arb_poly_mullow_block(res->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
    }

    _arb_poly_set_length(res, zlen);
    _arb_poly_normalise(res);
}

