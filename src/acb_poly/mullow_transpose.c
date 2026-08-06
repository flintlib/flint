/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void
_acb_poly_mulmid_transpose(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong nlo, slong nhi, slong prec)
{
    arb_ptr a, b, c, d, e, f, w;
    arb_ptr t;
    slong i;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    slong nlo2 = (len1 + len2 - 1) - nlo;

    if (len1 > nlo2)
    {
        slong trunc = len1 - nlo2;
        poly1 += trunc;
        len1 -= trunc;
        nlo -= trunc;
        nhi -= trunc;
    }

    if (len2 > nlo2)
    {
        slong trunc = len2 - nlo2;
        poly2 += trunc;
        len2 -= trunc;
        nlo -= trunc;
        nhi -= trunc;
    }

    w = flint_malloc(sizeof(arb_struct) * (2 * (len1 + len2 + (nhi - nlo))));
    a = w;
    b = a + len1;
    c = b + len1;
    d = c + len2;
    e = d + len2;
    f = e + (nhi - nlo);

    /* (e+fi) = (a+bi)(c+di) = (ac - bd) + (ad + bc)i */
    t = _arb_vec_init(nhi - nlo);

    for (i = 0; i < len1; i++)
    {
        a[i] = *acb_realref(poly1 + i);
        b[i] = *acb_imagref(poly1 + i);
    }

    for (i = 0; i < len2; i++)
    {
        c[i] = *acb_realref(poly2 + i);
        d[i] = *acb_imagref(poly2 + i);
    }

    for (i = 0; i < nhi - nlo; i++)
    {
        e[i] = *acb_realref(res + i);
        f[i] = *acb_imagref(res + i);
    }

    _arb_poly_mulmid(e, a, len1, c, len2, nlo, nhi, prec);
    _arb_poly_mulmid(t, b, len1, d, len2, nlo, nhi, prec);
    _arb_vec_sub(e, e, t, nhi - nlo, prec);

    _arb_poly_mulmid(f, a, len1, d, len2, nlo, nhi, prec);
    /* squaring */
    if (poly1 == poly2 && len1 == len2)
    {
        _arb_vec_scalar_mul_2exp_si(f, f, nhi - nlo, 1);
    }
    else
    {
        _arb_poly_mulmid(t, b, len1, c, len2, nlo, nhi, prec);
        _arb_vec_add(f, f, t, nhi - nlo, prec);
    }

    for (i = 0; i < nhi - nlo; i++)
    {
        *acb_realref(res + i) = e[i];
        *acb_imagref(res + i) = f[i];
    }

    _arb_vec_clear(t, nhi - nlo);
    flint_free(w);
}

void
_acb_poly_mullow_transpose(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec)
{
    _acb_poly_mulmid_transpose(res, poly1, len1, poly2, len2, 0, n, prec);
}

void
acb_poly_mulmid_transpose(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong nlo, slong nhi, slong prec)
{
    slong xlen, ylen, zlen;

    xlen = poly1->length;
    ylen = poly2->length;

    if (xlen == 0 || ylen == 0 || nlo >= FLINT_MIN(nhi, xlen + ylen - 1))
    {
        acb_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, xlen + ylen - 1);
    zlen = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        acb_poly_t tmp;
        acb_poly_init2(tmp, zlen);
        _acb_poly_mulmid_transpose(tmp->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, nlo, nhi, prec);
        acb_poly_swap(res, tmp);
        acb_poly_clear(tmp);
    }
    else
    {
        acb_poly_fit_length(res, zlen);
        _acb_poly_mulmid_transpose(res->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, nlo, nhi, prec);
    }

    _acb_poly_set_length(res, zlen);
    _acb_poly_normalise(res);
}

void
acb_poly_mullow_transpose(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                slong n, slong prec)
{
    slong len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    n = FLINT_MIN((len1 + len2 - 1), n);
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (res == poly1 || res == poly2)
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_mullow_transpose(t->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
    else
    {
        acb_poly_fit_length(res, n);
        _acb_poly_mullow_transpose(res->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}
