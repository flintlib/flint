/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_extras.h"
#include "fmpzi.h"
#include "arb_poly.h"
#include "acb_poly.h"

/*
    Multiplication of well-scaled complex polynomials through a single
    Gaussian integer polynomial product.

    After the substitution x -> 2^(-c) x (c chosen as in arb_poly's
    block multiplication, from the first and last nonzero midpoints),
    the real and imaginary midpoints of each operand are written as
    Gaussian integers on one common exponent; this is exact, and is
    done only when the height of the operand (top bit minus common
    bottom bit) is at most ALPHA prec + BETA, the bound arb_poly uses
    for one block of a real polynomial. The midpoints of the product are
    then one _fmpzi_poly_mulmid call, rounded to prec, against four real
    products for the transposed algorithm (which blocks each part on its
    own and so also handles poorly scaled input, to which we leave it).

    The radii are bounded separately, component by component: for
    (ar + i ai)(br + i bi) with radii rar, rai, rbr, rbi,

      rad re <= |ar| rbr + rar (|br| + rbr) + |ai| rbi + rai (|bi| + rbi),
      rad im <= |ar| rbi + rar (|bi| + rbi) + |ai| rbr + rai (|br| + rbr),

    summed over the convolution: up to eight products of nonnegative
    polynomials (four for squares), each skipped when a factor vanishes,
    computed by the radius machinery of arb_poly's block multiplication
    with the same scaling.
*/

#define ALPHA 3.0
#define BETA 512

static fmpzi_struct *
_fmpzi_vec_init(slong len)
{
    return (fmpzi_struct *) flint_calloc(len, sizeof(fmpzi_struct));
}

static void
_fmpzi_vec_clear(fmpzi_struct * v, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        fmpzi_clear(v + i);
    flint_free(v);
}

static int
_fmpzi_vec_is_zero(const fmpzi_struct * v, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        if (!fmpzi_is_zero(v + i))
            return 0;
    return 1;
}

/* exponent of the larger part of the midpoint of x; 0 if it is zero */
static int
_acb_mid_top(fmpz_t e, const acb_t x)
{
    const arf_struct * a = arb_midref(acb_realref(x));
    const arf_struct * b = arb_midref(acb_imagref(x));

    if (arf_is_zero(a) && arf_is_zero(b))
        return 0;

    if (arf_is_zero(b) || (!arf_is_zero(a) && fmpz_cmp(ARF_EXPREF(a), ARF_EXPREF(b)) >= 0))
        fmpz_set(e, ARF_EXPREF(a));
    else
        fmpz_set(e, ARF_EXPREF(b));

    return 1;
}

static int
_acb_mid_is_zero(const acb_t x)
{
    return arf_is_zero(arb_midref(acb_realref(x))) && arf_is_zero(arb_midref(acb_imagref(x)));
}

/* as _arb_poly_get_scale, with the larger part of each midpoint */
static void
_acb_poly_get_scale(fmpz_t scale, acb_srcptr x, slong xlen,
                                  acb_srcptr y, slong ylen)
{
    slong xa, xb, ya, yb, den;
    fmpz_t e;

    fmpz_zero(scale);

    xa = 0;
    xb = xlen - 1;
    while (xa < xlen && _acb_mid_is_zero(x + xa)) xa++;
    while (xb > xa && _acb_mid_is_zero(x + xb)) xb--;

    ya = 0;
    yb = ylen - 1;
    while (ya < ylen && _acb_mid_is_zero(y + ya)) ya++;
    while (yb > ya && _acb_mid_is_zero(y + yb)) yb--;

    if (xa <= xb && ya <= yb && (xa < xb || ya < yb))
    {
        fmpz_init(e);

        _acb_mid_top(e, x + xb);
        fmpz_add(scale, scale, e);
        _acb_mid_top(e, x + xa);
        fmpz_sub(scale, scale, e);
        _acb_mid_top(e, y + yb);
        fmpz_add(scale, scale, e);
        _acb_mid_top(e, y + ya);
        fmpz_sub(scale, scale, e);

        den = (xb - xa) + (yb - ya);

        /* scale = floor(scale / den + 1/2) */
        fmpz_mul_2exp(scale, scale, 1);
        fmpz_add_ui(scale, scale, den);
        fmpz_fdiv_q_ui(scale, scale, 2 * den);

        fmpz_clear(e);
    }
}

/*
    Writes the midpoints of (x, len), divided by 2^(scale i), as
    z_i 2^bot with Gaussian integers z_i on one exponent bot. Returns 0
    (leaving z undefined) if that takes more than maxheight bits.
*/
static int
_acb_vec_get_fmpzi_2exp(fmpzi_struct * z, fmpz_t bot, acb_srcptr x,
                        slong len, const fmpz_t scale, slong maxheight)
{
    fmpz_t top, t, b, si;
    slong i, j, s;
    int any = 0, ok = 1;

    fmpz_init(top);
    fmpz_init(t);
    fmpz_init(b);
    fmpz_init(si);

    /* the top and bottom bits */
    for (i = 0; i < len; i++)
    {
        fmpz_mul_si(si, scale, i);

        for (j = 0; j < 2; j++)
        {
            const arf_struct * m = arb_midref(j ? acb_imagref(x + i) : acb_realref(x + i));

            if (arf_is_zero(m))
                continue;

            fmpz_sub(t, ARF_EXPREF(m), si);
            fmpz_sub_ui(b, t, arf_bits(m));

            if (!any)
            {
                fmpz_swap(top, t);
                fmpz_swap(bot, b);
                any = 1;
            }
            else
            {
                if (fmpz_cmp(t, top) > 0)
                    fmpz_swap(top, t);
                if (fmpz_cmp(b, bot) < 0)
                    fmpz_swap(bot, b);
            }
        }
    }

    if (!any)
    {
        fmpz_zero(bot);
        for (i = 0; i < len; i++)
            fmpzi_zero(z + i);
        goto cleanup;
    }

    fmpz_sub(t, top, bot);
    if (fmpz_cmp_si(t, maxheight) > 0)
    {
        ok = 0;
        goto cleanup;
    }

    for (i = 0; i < len; i++)
    {
        fmpz_mul_si(si, scale, i);

        for (j = 0; j < 2; j++)
        {
            const arf_struct * m = arb_midref(j ? acb_imagref(x + i) : acb_realref(x + i));
            fmpz * c = j ? fmpzi_imagref(z + i) : fmpzi_realref(z + i);

            if (arf_is_zero(m))
            {
                fmpz_zero(c);
                continue;
            }

            /* m = c 2^b, and c 2^(b - scale i) = (c 2^s) 2^bot */
            arf_get_fmpz_2exp(c, b, m);
            fmpz_sub(t, b, si);
            s = _fmpz_sub_small(t, bot);
            fmpz_mul_2exp(c, c, s);
        }
    }

cleanup:
    fmpz_clear(top);
    fmpz_clear(t);
    fmpz_clear(b);
    fmpz_clear(si);
    return ok;
}

int
_acb_poly_mulmid_fmpzi(acb_ptr res, acb_srcptr x, slong xlen,
    acb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)
{
    fmpzi_struct * xz, * yz, * zz;
    mag_ptr xmr, xmi, xrr, xri, ypr, ypi, yrr, yri, rre, rim;
    fmpz_t scale, xbot, ybot, e;
    slong i, w, w_all, maxheight;
    int squaring, ok, xzero, yzero;
    double h;

    xlen = FLINT_MIN(xlen, nhi);
    ylen = FLINT_MIN(ylen, nhi);

    /* the window may extend beyond the product */
    w_all = nhi - nlo;
    if (nlo >= xlen + ylen - 1)
    {
        _acb_vec_zero(res, w_all);
        return 1;
    }
    nhi = FLINT_MIN(nhi, xlen + ylen - 1);

    /* high truncation of the inputs */
    {
        slong nlo2 = (xlen + ylen - 1) - nlo;

        if (xlen > nlo2)
        {
            slong trunc = xlen - nlo2;
            x += trunc;
            xlen -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (ylen > nlo2)
        {
            slong trunc = ylen - nlo2;
            y += trunc;
            ylen -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

    squaring = (x == y && xlen == ylen);

    if (!_acb_vec_is_finite(x, xlen) || (!squaring && !_acb_vec_is_finite(y, ylen)))
        return 0;

    h = ALPHA * (double) prec + BETA;
    maxheight = (h > (double) (WORD_MAX / 4)) ? WORD_MAX / 4 : (slong) h;

    w = nhi - nlo;

    fmpz_init(scale);
    fmpz_init(xbot);
    fmpz_init(ybot);
    fmpz_init(e);
    xz = _fmpzi_vec_init(xlen);
    yz = squaring ? xz : _fmpzi_vec_init(ylen);
    zz = NULL;

    _acb_poly_get_scale(scale, x, xlen, y, ylen);

    ok = _acb_vec_get_fmpzi_2exp(xz, xbot, x, xlen, scale, maxheight);
    if (ok && !squaring)
        ok = _acb_vec_get_fmpzi_2exp(yz, ybot, y, ylen, scale, maxheight);
    else
        fmpz_set(ybot, xbot);

    if (!ok)
        goto cleanup;

    /* the midpoints */
    xzero = _fmpzi_vec_is_zero(xz, xlen);
    yzero = squaring ? xzero : _fmpzi_vec_is_zero(yz, ylen);
    zz = _fmpzi_vec_init(w);
    if (!xzero && !yzero)
        _fmpzi_poly_mulmid(zz, xz, xlen, yz, ylen, nlo, nhi);

    /* the radii, before res is written (it may alias an input) */
    rre = _mag_vec_init(2 * w);
    rim = rre + w;

    xmr = _mag_vec_init(4 * xlen);
    xmi = xmr + xlen;
    xrr = xmi + xlen;
    xri = xrr + xlen;

    for (i = 0; i < xlen; i++)
    {
        mag_set(xrr + i, arb_radref(acb_realref(x + i)));
        mag_set(xri + i, arb_radref(acb_imagref(x + i)));
    }

    if (squaring)
    {
        /* rad re <= rar (2 |ar| + rar) + rai (2 |ai| + rai),
           rad im <= rai (2 |ar| + rar) + rar (2 |ai| + rai) */
        for (i = 0; i < xlen; i++)
        {
            arf_get_mag(xmr + i, arb_midref(acb_realref(x + i)));
            mag_mul_2exp_si(xmr + i, xmr + i, 1);
            mag_add(xmr + i, xmr + i, xrr + i);
            arf_get_mag(xmi + i, arb_midref(acb_imagref(x + i)));
            mag_mul_2exp_si(xmi + i, xmi + i, 1);
            mag_add(xmi + i, xmi + i, xri + i);
        }

        _arb_poly_addmulmid_mag(rre, xrr, xlen, xmr, xlen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rre, xri, xlen, xmi, xlen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rim, xri, xlen, xmr, xlen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rim, xrr, xlen, xmi, xlen, nlo, nhi, scale);
    }
    else
    {
        ypr = _mag_vec_init(4 * ylen);
        ypi = ypr + ylen;
        yrr = ypi + ylen;
        yri = yrr + ylen;

        for (i = 0; i < xlen; i++)
        {
            arf_get_mag(xmr + i, arb_midref(acb_realref(x + i)));
            arf_get_mag(xmi + i, arb_midref(acb_imagref(x + i)));
        }

        for (i = 0; i < ylen; i++)
        {
            arb_get_mag(ypr + i, acb_realref(y + i));
            arb_get_mag(ypi + i, acb_imagref(y + i));
            mag_set(yrr + i, arb_radref(acb_realref(y + i)));
            mag_set(yri + i, arb_radref(acb_imagref(y + i)));
        }

        _arb_poly_addmulmid_mag(rre, xmr, xlen, yrr, ylen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rre, xrr, xlen, ypr, ylen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rre, xmi, xlen, yri, ylen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rre, xri, xlen, ypi, ylen, nlo, nhi, scale);

        _arb_poly_addmulmid_mag(rim, xmr, xlen, yri, ylen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rim, xrr, xlen, ypi, ylen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rim, xmi, xlen, yrr, ylen, nlo, nhi, scale);
        _arb_poly_addmulmid_mag(rim, xri, xlen, ypr, ylen, nlo, nhi, scale);

        _mag_vec_clear(ypr, 4 * ylen);
    }

    _mag_vec_clear(xmr, 4 * xlen);

    /* res_k = zz_k 2^(xbot + ybot + scale k), rounded, plus the radii */
    fmpz_add(e, xbot, ybot);
    fmpz_addmul_ui(e, scale, nlo);

    for (i = 0; i < w; i++)
    {
        arb_set_round_fmpz_2exp(acb_realref(res + i), fmpzi_realref(zz + i), e, prec);
        arb_set_round_fmpz_2exp(acb_imagref(res + i), fmpzi_imagref(zz + i), e, prec);
        arb_add_error_mag(acb_realref(res + i), rre + i);
        arb_add_error_mag(acb_imagref(res + i), rim + i);
        fmpz_add(e, e, scale);
    }

    _mag_vec_clear(rre, 2 * w);
    _acb_vec_zero(res + w, w_all - w);

cleanup:
    _fmpzi_vec_clear(xz, xlen);
    if (!squaring)
        _fmpzi_vec_clear(yz, ylen);
    if (zz != NULL)
        _fmpzi_vec_clear(zz, w);
    fmpz_clear(scale);
    fmpz_clear(xbot);
    fmpz_clear(ybot);
    fmpz_clear(e);

    return ok;
}

/*
    Whether to use _acb_poly_mulmid_fmpzi rather than the classical
    algorithm. Measured (x86-64, random well-scaled coefficients with
    radii, products, low products and squares), the crossover is at a
    nearly constant number q of coefficient products in the window per
    input coefficient for each precision, about 20 at 64 bits, 10 at
    256 bits, 4.5 at 1024 bits and 1 at 4096 bits for products, i.e.
    q prec^0.57 >= 240, and about 1.4 times that for squares, whose
    classical algorithm halves the products.
*/
int
_acb_poly_mulmid_use_fmpzi(slong len1, slong len2, slong nlo, slong nhi,
                           slong prec, int squaring)
{
    double terms = 0, q;
    slong k;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);
    nhi = FLINT_MIN(nhi, len1 + len2 - 1);

    if (len1 <= 1 || len2 <= 1 || nlo >= nhi)
        return 0;

    for (k = nlo; k < nhi; k++)
        terms += FLINT_MIN(k, len1 - 1) - FLINT_MAX(0, k - len2 + 1) + 1;

    q = terms / (double) (len1 + len2);

    if (q < 1.0)
        return 0;

    return q * pow((double) FLINT_MAX(prec, 2), 0.57) >= (squaring ? 330.0 : 240.0);
}
