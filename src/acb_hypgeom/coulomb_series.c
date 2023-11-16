/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_hypgeom_coulomb_series(acb_ptr F, acb_ptr G,
    acb_ptr Hpos, acb_ptr Hneg, const acb_t l, const acb_t eta,
    acb_srcptr z, slong zlen, slong len, slong prec)
{
    acb_ptr t, v;

    if (len == 0)
        return;

    zlen = FLINT_MIN(zlen, len);

    if (zlen == 1)
    {
        acb_hypgeom_coulomb(F, G, Hpos, Hneg, l, eta, z, prec);
        if (F != NULL) _acb_vec_zero(F + 1, len - 1);
        if (G != NULL) _acb_vec_zero(G + 1, len - 1);
        if (Hpos != NULL) _acb_vec_zero(Hpos + 1, len - 1);
        if (Hneg != NULL) _acb_vec_zero(Hneg + 1, len - 1);
        return;
    }

    t = _acb_vec_init(len);
    v = _acb_vec_init(zlen);

    /* copy nonconstant part first to allow aliasing */
    acb_zero(v);
    _acb_vec_set(v + 1, z + 1, zlen - 1);

    acb_hypgeom_coulomb_jet(F, G, Hpos, Hneg, l, eta, z, len, prec);

    if (F != NULL)
    {
        _acb_vec_set(t, F, len);
        _acb_poly_compose_series(F, t, len, v, zlen, len, prec);
    }

    if (G != NULL)
    {
        _acb_vec_set(t, G, len);
        _acb_poly_compose_series(G, t, len, v, zlen, len, prec);
    }

    if (Hpos != NULL)
    {
        _acb_vec_set(t, Hpos, len);
        _acb_poly_compose_series(Hpos, t, len, v, zlen, len, prec);
    }

    if (Hneg != NULL)
    {
        _acb_vec_set(t, Hneg, len);
        _acb_poly_compose_series(Hneg, t, len, v, zlen, len, prec);
    }

    _acb_vec_clear(t, len);
    _acb_vec_clear(v, zlen);
}

void
acb_hypgeom_coulomb_series(acb_poly_t F, acb_poly_t G,
    acb_poly_t Hpos, acb_poly_t Hneg, const acb_t l, const acb_t eta,
        const acb_poly_t z, slong len, slong prec)
{
    acb_srcptr zptr;
    slong zlen;
    acb_t t;

    if (len == 0)
    {
        if (F != NULL) acb_poly_zero(F);
        if (G != NULL) acb_poly_zero(G);
        if (Hpos != NULL) acb_poly_zero(Hpos);
        if (Hneg != NULL) acb_poly_zero(Hneg);
        return;
    }

    zlen = z->length;
    if (zlen <= 1)
        len = 1;

    if (F != NULL) acb_poly_fit_length(F, len);
    if (G != NULL) acb_poly_fit_length(G, len);
    if (Hpos != NULL) acb_poly_fit_length(Hpos, len);
    if (Hneg != NULL) acb_poly_fit_length(Hneg, len);

    if (zlen == 0)
    {
        acb_init(t);
        zptr = t;
        zlen = 1;
    }
    else
    {
        zptr = z->coeffs;
    }

    _acb_hypgeom_coulomb_series(
        F ? F->coeffs : NULL,
        G ? G->coeffs : NULL,
        Hpos ? Hpos->coeffs : NULL,
        Hneg ? Hneg->coeffs : NULL,
        l, eta,
        zptr, zlen, len, prec);

    if (F != NULL) _acb_poly_set_length(F, len);
    if (G != NULL) _acb_poly_set_length(G, len);
    if (Hpos != NULL) _acb_poly_set_length(Hpos, len);
    if (Hneg != NULL) _acb_poly_set_length(Hneg, len);

    if (F != NULL) _acb_poly_normalise(F);
    if (G != NULL) _acb_poly_normalise(G);
    if (Hpos != NULL) _acb_poly_normalise(Hpos);
    if (Hneg != NULL) _acb_poly_normalise(Hneg);
}

