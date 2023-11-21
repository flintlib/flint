/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_hypgeom_airy_series(acb_ptr ai, acb_ptr ai_prime,
    acb_ptr bi, acb_ptr bi_prime, acb_srcptr z, slong zlen, slong len, slong prec)
{
    acb_ptr t, u, v;
    slong tlen = len + ((ai_prime != NULL) || (bi_prime != NULL));

    zlen = FLINT_MIN(zlen, len);

    if (zlen <= 0)
        return;

    if (zlen == 1)
    {
        acb_hypgeom_airy(ai, ai_prime, bi, bi_prime, z, prec);
        return;
    }

    t = _acb_vec_init(tlen);
    u = _acb_vec_init(tlen);
    v = _acb_vec_init(len);

    acb_hypgeom_airy_jet((ai || ai_prime) ? t : NULL,
                         (bi || bi_prime) ? u : NULL, z, tlen, prec);

    /* compose with nonconstant part */
    acb_zero(v);
    _acb_vec_set(v + 1, z + 1, zlen - 1);

    if (ai != NULL) _acb_poly_compose_series(ai, t, len, v, zlen, len, prec);
    if (bi != NULL) _acb_poly_compose_series(bi, u, len, v, zlen, len, prec);

    /* todo: use chain rule to avoid compositions for derivatives? */
    if (ai_prime != NULL)
    {
        _acb_poly_derivative(t, t, len + 1, prec);
        _acb_poly_compose_series(ai_prime, t, len, v, zlen, len, prec);
    }

    if (bi_prime != NULL)
    {
        _acb_poly_derivative(u, u, len + 1, prec);
        _acb_poly_compose_series(bi_prime, u, len, v, zlen, len, prec);
    }

    _acb_vec_clear(t, tlen);
    _acb_vec_clear(u, tlen);
    _acb_vec_clear(v, len);
}

void
acb_hypgeom_airy_series(acb_poly_t ai, acb_poly_t ai_prime,
    acb_poly_t bi, acb_poly_t bi_prime, const acb_poly_t z, slong len, slong prec)
{
    if (len == 0)
    {
        if (ai != NULL) acb_poly_zero(ai);
        if (ai_prime != NULL) acb_poly_zero(ai_prime);
        if (bi != NULL) acb_poly_zero(bi);
        if (bi_prime != NULL) acb_poly_zero(bi_prime);
        return;
    }

    if (z->length <= 1)
        len = 1;

    if (ai != NULL) acb_poly_fit_length(ai, len);
    if (ai_prime != NULL) acb_poly_fit_length(ai_prime, len);
    if (bi != NULL) acb_poly_fit_length(bi, len);
    if (bi_prime != NULL) acb_poly_fit_length(bi_prime, len);

    if (z->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_hypgeom_airy_series(
            ai ? ai->coeffs : NULL, ai_prime ? ai_prime->coeffs : NULL,
            bi ? bi->coeffs : NULL, bi_prime ? bi_prime->coeffs : NULL,
            t, 1, len, prec);
        acb_clear(t);
    }
    else
    {
        _acb_hypgeom_airy_series(
            ai ? ai->coeffs : NULL, ai_prime ? ai_prime->coeffs : NULL,
            bi ? bi->coeffs : NULL, bi_prime ? bi_prime->coeffs : NULL,
            z->coeffs, z->length, len, prec);
    }

    if (ai != NULL)
    {
        _acb_poly_set_length(ai, len);
        _acb_poly_normalise(ai);
    }

    if (ai_prime != NULL)
    {
        _acb_poly_set_length(ai_prime, len);
        _acb_poly_normalise(ai_prime);
    }

    if (bi != NULL)
    {
        _acb_poly_set_length(bi, len);
        _acb_poly_normalise(bi);
    }

    if (bi_prime != NULL)
    {
        _acb_poly_set_length(bi_prime, len);
        _acb_poly_normalise(bi_prime);
    }
}

