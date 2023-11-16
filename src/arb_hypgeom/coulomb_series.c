/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"

void
_arb_hypgeom_coulomb_series(arb_ptr F, arb_ptr G,
    const arb_t l, const arb_t eta,
    arb_srcptr z, slong zlen, slong len, slong prec)
{
    arb_ptr t, v;

    if (len == 0)
        return;

    zlen = FLINT_MIN(zlen, len);

    if (zlen == 1)
    {
        arb_hypgeom_coulomb(F, G, l, eta, z, prec);
        if (F != NULL) _arb_vec_zero(F + 1, len - 1);
        if (G != NULL) _arb_vec_zero(G + 1, len - 1);
        return;
    }

    t = _arb_vec_init(len);
    v = _arb_vec_init(zlen);

    /* copy nonconstant part first to allow aliasing */
    arb_zero(v);
    _arb_vec_set(v + 1, z + 1, zlen - 1);

    arb_hypgeom_coulomb_jet(F, G, l, eta, z, len, prec);

    if (F != NULL)
    {
        _arb_vec_set(t, F, len);
        _arb_poly_compose_series(F, t, len, v, zlen, len, prec);
    }

    if (G != NULL)
    {
        _arb_vec_set(t, G, len);
        _arb_poly_compose_series(G, t, len, v, zlen, len, prec);
    }

    _arb_vec_clear(t, len);
    _arb_vec_clear(v, zlen);
}

void
arb_hypgeom_coulomb_series(arb_poly_t F, arb_poly_t G,
        const arb_t l, const arb_t eta,
        const arb_poly_t z, slong len, slong prec)
{
    arb_srcptr zptr;
    slong zlen;
    arb_t t;

    if (len == 0)
    {
        if (F != NULL) arb_poly_zero(F);
        if (G != NULL) arb_poly_zero(G);
        return;
    }

    zlen = z->length;
    if (zlen <= 1)
        len = 1;

    if (F != NULL) arb_poly_fit_length(F, len);
    if (G != NULL) arb_poly_fit_length(G, len);

    if (zlen == 0)
    {
        arb_init(t);
        zptr = t;
        zlen = 1;
    }
    else
    {
        zptr = z->coeffs;
    }

    _arb_hypgeom_coulomb_series(
        F ? F->coeffs : NULL,
        G ? G->coeffs : NULL,
        l, eta,
        zptr, zlen, len, prec);

    if (F != NULL) _arb_poly_set_length(F, len);
    if (G != NULL) _arb_poly_set_length(G, len);

    if (F != NULL) _arb_poly_normalise(F);
    if (G != NULL) _arb_poly_normalise(G);
}

